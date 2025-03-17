#!/usr/bin/env python3
import os
import sys
import logging
import argparse
import concurrent.futures
import subprocess
import multiprocessing
from pathlib import Path
import shutil
import duckdb
import pandas as pd
from typing import List, Optional, Union
import pysam
import psutil
import shutil

# Add barcodes to BAM files
def add_barcode(bam_file, output_bam_bc_dir, threads, logger):
    # Probably not the most effiencent function
    """
    Add cell barcode to BAM file from filename.   
    Args:
        bam_file: Path to the BAM file
        output_bam_bc_dir: Output directory for BAM files that are extended with a barcode
        threads: Number of threads to use for indexing
        logger: Logger instance
    
    Returns:
        bool: Success status
    """
    try:
        # Extract filename without extension
        base_name = os.path.basename(bam_file)[:-4]
        cell_barcode = base_name
        output_path = os.path.join(output_bam_bc_dir, f"{base_name}.cb.bam")
        
        logger.info(f"Processing {bam_file} with barcode {cell_barcode}...")
        
        # Create an index file if it doesn't exist
        index_file = os.path.join(output_bam_bc_dir, f"{base_name}.bai")
        if not os.path.exists(index_file):
            logger.info(f"Indexing {bam_file}...")


            samtools_path = shutil.which("samtools")
            result = subprocess.run(
                [samtools_path, "index", "--threads", str(min(4, threads)), bam_file, index_file],
                check=True,
                stderr=subprocess.PIPE
            )
            if result.returncode != 0:
                logger.error(f"Failed to index {bam_file}: {result.stderr.decode()}")
                return False
        
        # Add cell barcode using pysam      
        logger.info(f"Adding barcode to {bam_file}...")
        infile = pysam.AlignmentFile(bam_file, "rb")
        outfile = pysam.AlignmentFile(output_path, "wb", header=infile.header)
        
        for read in infile:
            read.set_tag("CB", cell_barcode, "Z", replace=True)
            outfile.write(read)
            
        infile.close()
        outfile.close()
        
        # Remove the index file as we'll need to create a new one with the added barcodes for
        # the fragment file creation
        if os.path.exists(index_file):
            os.remove(index_file)
            
        logger.info(f"Successfully processed {bam_file}")
        return True
    
    except Exception as e:
        logger.error(f"Error processing {bam_file}: {str(e)}")
        return False

# Function to create fragment files from BAM files
def create_fragments(bam_file, output_bed_dir, max_procs_per_file=22, max_tn5_distance=1000, min_mapping_quality=30, logger=None):
    """
    Create fragment file from BAM file with cell barcodes. This is just a parallel file processing wraper around sinto.
    Args:
        bam_file: Path to the BAM file with barcodes
        output_bed_dir: Output directory for BED files
        max_procs_per_file: Maximum number of processes to use, one process requieres approcimately 1-2GB memory
        max_tn5_distance: Maximum TN5 distance
        min_mapping_quality: Minimum mapping quality
        logger: Logger instance
    
    Returns:
        bool: Success status
    """
    try:
        # Extract filename without extension, handle both .bam and .cb.bam extensions
        base_name = os.path.basename(bam_file)
        if base_name.endswith('.cb.bam'):
            base_name = base_name[:-7]
        elif base_name.endswith('.bam'):
            base_name = base_name[:-4]
            
        fragments_file = os.path.join(output_bed_dir, f"{base_name}.bed")
        
        logger.info(f"Creating fragment file for {bam_file}...")
        
        # Create index if it doesn't exist
        if not os.path.exists(f"{bam_file}.bai"):
            logger.info(f"Indexing {bam_file}...")
            result = subprocess.run(
                ["samtools", "index", "--threads", str(max_procs_per_file), bam_file],
                check=True,
                stderr=subprocess.PIPE
            )
            if result.returncode != 0:
                logger.error(f"Failed to index {bam_file}: {result.stderr.decode()}")
                return False
        
        # Create fragment file using sinto
        logger.info(f"Running sinto fragments for {bam_file} with {max_procs_per_file} threads...")
        sinto_path = shutil.which("sinto")
        result = subprocess.run([
            sinto_path, "fragments", 
            "-b", bam_file, 
            "-f", fragments_file, 
            "-p", str(max_procs_per_file),
            "--max_distance", str(max_tn5_distance),
            "--min_mapq", str(min_mapping_quality)
        ], check=True, stderr=subprocess.PIPE)
        
        if result.returncode != 0:
            logger.error(f"Failed to create fragment file for {bam_file}: {result.stderr.decode()}")
            return False
        
        # Remove the index file to save space
        if os.path.exists(f"{bam_file}.bai"):
            os.remove(f"{bam_file}.bai")
            
        logger.info(f"Successfully created fragment file for {bam_file}")
        return True
    
    except Exception as e:
        logger.error(f"Error creating fragment file for {bam_file}: {str(e)}")
        return False

# Function to process fragments with DuckDB
def process_fragments_with_duckdb(output_bed_dir, db_path, min_size=0, max_size=1000, 
                                  memory_limit='8GB', count_table_path=None, logger=None):
    """
    Process fragment files using DuckDB for efficient handling of large datasets.
    
    Args:
        output_bed_dir: Directory containing BED files
        db_path: Path to DuckDB database file
        min_size: Minimum fragment size
        max_size: Maximum fragment size
        memory_limit: Memory limit for DuckDB
        count_table_path: Path to save the count table
        logger: Logger instance
    
    Returns:
        bool: Success status
    """
    try:
        if logger:
            logger.info(f"Processing fragments with DuckDB to {db_path}...")
        
        # Get all fragment files
        fragment_files = list(map(str, Path(output_bed_dir).glob("*.bed")))
        
        if not fragment_files:
            if logger:
                logger.error("No fragment files found in the output directory")
            return False
        
        # Process fragments and insert into persistent DuckDB database
        insert_bed_to_duckdb(
            fragment_files=fragment_files,
            db_path=db_path,
            summarize=True, # Create the count_table that is used in the PeakQC scoring function
            min_size=min_size,
            max_size=max_size,
            memory_limit=memory_limit,
            count_table_path=count_table_path
        )
               
        if logger:
            logger.info(f"Successfully processed and concatenated fragment files using DuckDB")
        return True
    
    except Exception as e:
        if logger:
            logger.error(f"Error processing fragments with DuckDB: {str(e)}")
        return False

# Function to insert fragment BED files into DuckDB
def insert_bed_to_duckdb(
  fragment_files: List[str],
  db_path: str,
  summarize: bool = False,
  min_size: int = 0,
  max_size: int = 1000, 
  memory_limit: str = '8GB',
  count_table_path: Optional[str | Path] = None) -> pd.DataFrame:
    """
    Insert data from a fragment .bed file into the DuckDB database and optionally summarize the data. summarize = True will count the number of fragments, calculate the mean fragment length, and creates a fragment length distribution array for the respective barcode. The summarized data is stored inside the count_table table and the whole data is stored in the fragments table.
    
    Args:
        fragment_files (List[str]): List of paths to the .bed files.
        db_path (str): Path to the DuckDB database where data will be stored.
        summarize (bool): Whether to summarize the data while inserting it. Default is False.
        memory_limit: Memory limit for DuckDB.
        min_size: Minimum size threshold for fragments
        max_size: Maximum size threshold for fragments
        count_table_path: Optional path to save the count table as TSV.
    Returns:
        count_table: DataFrame with barcode statistics and distributions
    """
    con = duckdb.connect(db_path)
    con.execute(f"SET memory_limit='{memory_limit}'")
    
    con.execute("""
    CREATE TABLE IF NOT EXISTS fragments (
        chrom STRING,
        start INTEGER,
        "end" INTEGER,
        barcode STRING,
        size INTEGER,
        count INTEGER
    );
    """)
    
    if summarize:
      con.execute("""
      CREATE TABLE IF NOT EXISTS count_table (
          barcode STRING,
          insertsize_count INTEGER,
          mean_insertsize DOUBLE,
          dist VARCHAR
      );
      """)

    
    for fragment_file in fragment_files:
      con.execute("""
      INSERT INTO fragments
      SELECT
          column0 AS chrom,
          CAST(column1 AS INTEGER) AS start,
          CAST(column2 AS INTEGER) AS "end",
          column3 AS barcode,
          (CAST(column2 AS INTEGER) - CAST(column1 AS INTEGER)) AS size,
          CAST(column4 AS INTEGER) AS count
      FROM read_csv_auto(?, delim='\t', header=False);
      """, [fragment_file])

  
    if summarize:
        con.execute("""
        INSERT INTO count_table
        WITH fragment_sizes AS (
            SELECT barcode, size, count
            FROM fragments
            WHERE size BETWEEN ? AND ?
            AND barcode IN (SELECT DISTINCT barcode FROM read_csv_auto(?, delim='\t', header=False))
        ),
        barcode_stats AS (
            SELECT barcode,
                   SUM(count) AS insertsize_count,
                   AVG(size) AS mean_insertsize
            FROM fragment_sizes
            GROUP BY barcode
        ),
        size_counts AS (
            SELECT 
                barcode,
                size,
                COUNT(*) AS size_count
            FROM fragment_sizes
            GROUP BY barcode, size
        ),
        size_arrays AS (
            SELECT 
                sc.barcode,
                ARRAY_AGG(COALESCE(sc.size_count, 0) ORDER BY sc.size) AS dist
            FROM size_counts sc
            GROUP BY sc.barcode
        )
        SELECT 
            bs.*, 
            sa.dist
        FROM barcode_stats bs
        JOIN size_arrays sa ON bs.barcode = sa.barcode;
        """, [min_size, max_size, fragment_file])

    con.execute("""
    CREATE INDEX IF NOT EXISTS idx_barcode ON fragments (barcode);
    """)
    
    con.execute("""
    CREATE INDEX IF NOT EXISTS idx_barcode ON count_table (barcode);
    """)

    count_table = con.execute("""
    SELECT * FROM count_table;
    """).df()
    

    count_table.set_index('barcode', inplace=True)
    

    if count_table_path is not None:
        output_df = count_table.copy()
        
        #output_df['dist'] = output_df['dist'].apply(lambda x: ','.join(map(str, x)))
        
        output_df.to_csv(count_table_path, sep='\t', index=True)
    

    con.close()
    
    return count_table

def main():
    parser = argparse.ArgumentParser(description="Single Cell and Bulk ATAC-seq processing pipeline with DuckDB. This pipeline prepares BAM files for quality control with peakqc.add_fld_scoring")
    parser.add_argument("--input_bam_list", required=True, help="File containing list of input BAM files")
    parser.add_argument("--output_bed", required=True, help="Output directory for BED files")
    parser.add_argument("--output_bam_bc", help="Output directory for BAM files with barcodes (required if not skipping barcode adding step and if the files should not be overwritten)")
    parser.add_argument("--output_error", required=True, help="Output directory for error logs")
    parser.add_argument("--skip_barcode_step", action="store_true", help="Skip barcode addition step (if BAM files already have barcodes)")
    parser.add_argument("--max_tn5_distance", type=int, default=1000, help="Maximum TN5 distance")
    parser.add_argument("--min_mapping_quality", type=int, default=30, help="Minimum mapping quality")
    parser.add_argument("--barcode_workers", type=int, default=0, 
                        help="Number of parallel workers for barcode addition")
    parser.add_argument("--fragment_workers", type=int, default=0, 
                        help="Number of parallel workers for fragment creation")
    parser.add_argument("--threads_per_file", type=int, default=0, 
                        help="Number of threads per file for fragment creation")
    parser.add_argument("--duckdb_path", default="fragments.duckdb", 
                        help="Path to DuckDB database file")
    parser.add_argument("--memory_limit", default="8GB", 
                        help="Memory limit for DuckDB")
    parser.add_argument("--min_fragment_size", type=int, default=0, 
                        help="Minimum fragment size")
    parser.add_argument("--max_fragment_size", type=int, default=1000, 
                        help="Maximum fragment size")
    parser.add_argument("--count_table_path", default=None, 
                        help="Path to save count table with fragment files summarization. This file is required for the peakqc.add_fld_scoring function")
    args = parser.parse_args()
    
    import time
    start_time = time.time()

    
    # Validate arguments
    if not args.skip_barcode_step and not args.output_bam_bc:
        parser.error("--output_bam_bc is required when not using --skip_barcode_step")
    
    # Create output directories
    directories_to_create = [args.output_bed, args.output_error]
    if not args.skip_barcode_step:
        directories_to_create.append(args.output_bam_bc)
        
    for directory in directories_to_create:
        os.makedirs(directory, exist_ok=True)
    
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(os.path.join(args.output_error, "error.log")),
            logging.StreamHandler()
        ]
    )
    logger = logging.getLogger(__name__)
    
    # Determine available processors and memory
    available_procs = multiprocessing.cpu_count()
    
    
       
    # Read input BAM files
    with open(args.input_bam_list, 'r') as f:
        bam_files = [line.strip() for line in f if line.strip()]
        
    #####################################################################################################
    bam_files = bam_files[:10]
    #####################################################################################################
    
    logger.info(f"Total BAM files to process: {len(bam_files)}")
    
    if args.skip_barcode_step:
        logger.info("Skipping barcode addition step (assuming BAM files already have barcodes)")
        if args.output_bam_bc:
          bam_files_for_fragments = list(Path(args.output_bam_bc).glob("*.cb.bam"))
        else:
          bam_files_for_fragments = bam_files
    else:
        barcode_start_time = time.time()

        barcode_workers = args.barcode_workers if args.barcode_workers > 0 else min(available_procs, 20)
        threads_for_barcodes = max(1, min(4, available_procs // max(1, barcode_workers)))
        
        logger.info(f"Using {barcode_workers} parallel workers for barcode addition")
        logger.info(f"Using {threads_for_barcodes} threads per file for barcode operations")
        

        logger.info("Adding cell barcode tags...")
        
        barcode_workers = min(barcode_workers, len(bam_files))
        
        with concurrent.futures.ProcessPoolExecutor(max_workers=barcode_workers) as executor:
            futures = [
                executor.submit(
                    add_barcode, 
                    bam_file, 
                    args.output_bam_bc, 
                    threads_for_barcodes,
                    logger
                ) 
                for bam_file in bam_files
            ]
            
            for future in concurrent.futures.as_completed(futures):
                try:
                    future.result()
                except Exception as e:
                    logger.error(f"Error in barcode addition task: {str(e)}")
        
        # Get all bam files with added barcodes
        bam_files_for_fragments = list(Path(args.output_bam_bc).glob("*.cb.bam"))
        
        barcode_end_time = time.time()
        logger.info(f"Barcode addition completed in {barcode_end_time - barcode_start_time:.2f} seconds")

    

    memory_per_thread_gb = 2
    total_memory_gb = psutil.virtual_memory().total / (1024**3)


    if args.threads_per_file > 0:
        threads_per_file = args.threads_per_file
    else:
        threads_per_file = min(22, available_procs, max(1, int(total_memory_gb / memory_per_thread_gb)))
    
    if args.fragment_workers > 0:
        fragment_workers = args.fragment_workers
    else:

        max_parallel_by_memory = max(1, int(total_memory_gb / (threads_per_file * memory_per_thread_gb)))

        max_parallel_by_cpu = max(1, available_procs // threads_per_file)
        fragment_workers = min(max_parallel_by_memory, max_parallel_by_cpu)
    
    logger.info("Creating fragment files...")
    fragment_start_time = time.time()


    fragment_workers = min(fragment_workers, len(bam_files_for_fragments))
    
    logger.info(f"Using {fragment_workers} parallel workers for fragment creation")
    logger.info(f"Using {threads_per_file} threads per file for fragment creation")
    logger.info(f"Processing {fragment_workers} BAM files in parallel for fragment creation")
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=fragment_workers) as executor:
        futures = [
            executor.submit(
                create_fragments, 
                str(bam_file), 
                args.output_bed, 
                threads_per_file,
                args.max_tn5_distance, 
                args.min_mapping_quality, 
                logger
            ) 
            for bam_file in bam_files_for_fragments
        ]
        
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as e:
                logger.error(f"Error in fragment creation task: {str(e)}")
    
    fragment_end_time = time.time()
    logger.info(f"Fragment creation completed in {fragment_end_time - fragment_start_time:.2f} seconds")


    count_table_start_time = time.time()
    logger.info("Processing fragments with DuckDB...")
    process_fragments_with_duckdb(
        output_bed_dir=args.output_bed,
        db_path=args.duckdb_path,
        min_size=args.min_fragment_size,
        max_size=args.max_fragment_size,
        memory_limit=args.memory_limit,
        count_table_path=args.count_table_path,
        logger=logger
    )
    
    count_table_end_time = time.time()
    logger.info(f"Count table and DuckDB database creation completed in {count_table_end_time - count_table_start_time:.2f} seconds")
    
    logger.info("Pipeline completed successfully")
    end_time = time.time()
    logger.info(f"Total execution time: {end_time - start_time:.2f} seconds")

if __name__ == "__main__":
    main()