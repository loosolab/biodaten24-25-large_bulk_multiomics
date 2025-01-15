import duckdb
import pandas as pd
import numpy as np
from typing import Optional
from beartype import beartype
from pathlib import Path

@beartype
def insertsize_from_fragments(
    fragments: str | Path,
    min_size: int = 0,
    max_size: int = 1000,
    memory_limit: str = '4GB'
) -> pd.DataFrame:
    """
    Process fragment file and calculate size distributions per barcode.
    Uses DuckDB's native streaming capabilities to handle large files.
    
    Args:
        fragments: Path to input fragment file
        min_size: Minimum size threshold
        max_size: Maximum size threshold
        memory_limit: Memory limit for DuckDB
    
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
    con.close()
    
    return count_table