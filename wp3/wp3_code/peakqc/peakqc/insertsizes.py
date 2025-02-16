import duckdb
import pandas as pd
import numpy as np
from typing import List, Optional, Union
from beartype import beartype
from pathlib import Path


@beartype
def insertsize_from_fragments(
    fragments: str | Path,
    min_size: int = 0,
    max_size: int = 1000,
    memory_limit: str = '8GB',
    count_table_path: Optional[str | Path] = None
) -> pd.DataFrame:
    """
    Process fragment file and calculate size distributions per barcode.
    
    Args:
        fragments: Path to input fragment file
        min_size: Minimum size threshold
        max_size: Maximum size threshold
        memory_limit: Memory limit for DuckDB
        count_table_path: Optional path to save the count table as TSV

    Returns:
        count_table: DataFrame with barcode statistics and distributions
    """
    con = duckdb.connect()
    con.execute(f"SET memory_limit='{memory_limit}'")
    
    count_table = con.execute("""
        WITH fragment_sizes AS (
            SELECT 
                column3 AS barcode,
                (CAST(column2 AS INTEGER) - CAST(column1 AS INTEGER) - 9) AS size,
                CAST(column4 AS INTEGER) as count
            FROM read_csv_auto(?, delim='\t', header=False)
            WHERE (CAST(column2 AS INTEGER) - CAST(column1 AS INTEGER) - 9) BETWEEN ? AND ?
        ),
        barcode_stats AS (
            SELECT 
                barcode,
                SUM(count) as insertsize_count,
                AVG(size) as mean_insertsize
            FROM fragment_sizes
            GROUP BY barcode
        ),
        all_combinations AS (
            SELECT 
                b.barcode,
                s.size
            FROM (SELECT DISTINCT barcode FROM fragment_sizes) b
            CROSS JOIN generate_series(?, ?) s(size)
        ),
        size_counts AS (
            SELECT 
                barcode,
                size,
                COUNT(*) as fragment_count
            FROM fragment_sizes
            GROUP BY barcode, size
        ),
        size_arrays AS (
            SELECT 
                ac.barcode,
                ARRAY_AGG(COALESCE(sc.fragment_count, 0) ORDER BY ac.size) AS dist
            FROM all_combinations ac
            LEFT JOIN size_counts sc 
                ON ac.barcode = sc.barcode 
                AND ac.size = sc.size
            GROUP BY ac.barcode
        )
        SELECT 
            bs.*,
            sa.dist
        FROM barcode_stats bs
        JOIN size_arrays sa ON bs.barcode = sa.barcode
    """, [fragments, min_size, max_size, min_size, max_size]).df()
    
    count_table.set_index('barcode', inplace=True)


    if count_table_path is not None:
        output_df = count_table.copy()
        output_df['dist'] = output_df['dist'].apply(lambda x: ','.join(map(str, x)))
        output_df.to_csv(count_table_path, sep='\t', index = True)

    con.close()
    
    return count_table
    

@beartype
def insert_bed_to_duckdb(
  fragment_files: List[str],
  db_path: str,
  summarize: bool = False,
  memory_limit: str = '8GB',
  count_table_path: Optional[str | Path] = None) -> pd.DataFrame:
    """
    Insert data from a fragment .bed file into the DuckDB database and optionally summarize the data. summarize = True will count the number of fragments, calculate the mean fragment length, and creates a fragment length distribution array for the respective barcode. The summarized data is stored inside the count_table table and the whole data is stored in the fragments table.
    
    Args:
        fragment_files (List[str]): List of paths to the .bed files.
        db_path (str): Path to the DuckDB database where data will be stored.
        summarize (bool): Whether to summarize the data while inserting it. Default is False.
        memory_limit: Memory limit for DuckDB.
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
        end INTEGER,
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
          dist ARRAY<INTEGER>
      );
      """)

    
    for fragment_file in fragment_files:
      con.execute("""
      INSERT INTO fragments
      SELECT
          column1 AS chrom,
          CAST(column2 AS INTEGER) AS start,
          CAST(column3 AS INTEGER) AS end,
          column4 AS barcode,
          (CAST(column3 AS INTEGER) - CAST(column2 AS INTEGER)) AS size,
          CAST(column5 AS INTEGER) AS count
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
        
        output_df['dist'] = output_df['dist'].apply(lambda x: ','.join(map(str, x)))
        
        output_df.to_csv(count_table_path, sep='\t', index=True)
    

    con.close()
    
    return count_table
    

def insertsize_from_duckdb(
    db_path: str, 
    min_size: int = 0, 
    max_size: int = 1000, 
    memory_limit: str = '8GB', 
    count_table_path: Optional[str | Path] = None
) -> pd.DataFrame:
    """
    Summarizes all fragments in the DuckDB database.
    Calculates size distributions per barcode from the fragments table in the database.
    
    Args:
        db_path: Path to the DuckDB database
        min_size: Minimum size threshold for fragments
        max_size: Maximum size threshold for fragments
        memory_limit: Memory limit for DuckDB
        count_table_path: Optional path to save the summary table as TSV

    Returns:
        count_table: DataFrame with barcode statistics and distributions
    """

    con = duckdb.connect(db_path)
    con.execute(f"SET memory_limit='{memory_limit}'")


    count_table = con.execute("""
        WITH fragment_sizes AS (
            SELECT 
                barcode,
                size,
                count
            FROM fragments
            WHERE size BETWEEN ? AND ?
        ),
        barcode_stats AS (
            SELECT 
                barcode,
                SUM(count) AS insertsize_count,
                AVG(size) AS mean_insertsize
            FROM fragment_sizes
            GROUP BY barcode
        ),
        all_combinations AS (
            SELECT 
                b.barcode,
                s.size
            FROM (SELECT DISTINCT barcode FROM fragment_sizes) b
            CROSS JOIN generate_series(?, ?) s(size)
        ),
        size_counts AS (
            SELECT 
                barcode,
                size,
                COUNT(*) AS fragment_count
            FROM fragment_sizes
            GROUP BY barcode, size
        ),
        size_arrays AS (
            SELECT 
                ac.barcode,
                ARRAY_AGG(COALESCE(sc.fragment_count, 0) ORDER BY ac.size) AS dist
            FROM all_combinations ac
            LEFT JOIN size_counts sc 
                ON ac.barcode = sc.barcode 
                AND ac.size = sc.size
            GROUP BY ac.barcode
        )
        SELECT 
            bs.*,
            sa.dist
        FROM barcode_stats bs
        JOIN size_arrays sa ON bs.barcode = sa.barcode
    """, [min_size, max_size, min_size, max_size]).df()

    count_table.set_index('barcode', inplace=True)

    if count_table_path is not None:
        output_df = count_table.copy()
        output_df['dist'] = output_df['dist'].apply(lambda x: ','.join(map(str, x)))
        output_df.to_csv(count_table_path, sep='\t', index=True)

    con.close()

    return count_table


