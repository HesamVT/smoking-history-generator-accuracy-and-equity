"""
File: shg_aggr_results.py
Author: Hesam Mahmoudi
Affiliation: MGH Center for Health Technology Assessment,
    Mass General Hospital, Harvard Medical School
Created: 2026-03-25
Project Description: Aggregate cohort-specific smoking-history summary tables
    into manuscript-ready SVG and PNG outputs.
"""

from __future__ import annotations

from datetime import datetime, timezone
from html import escape
from html.parser import HTMLParser
from pathlib import Path
import re
from textwrap import wrap
from xml.sax.saxutils import escape as xml_escape
import zipfile

import matplotlib.pyplot as plt
from matplotlib import transforms
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = Path(__file__).resolve().parent

# ---------------------------------------------------------------------------
# Input files and output targets
# ---------------------------------------------------------------------------

MEDIAN_CSV_PATHS: dict[str, Path] = {
    "SCCS": ROOT / "SCCS data" / "SCCS_median_tables_long.csv",
    "MECS": ROOT / "MECS data" / "MECS_median_tables_long.csv",
    "BRFSS": ROOT / "BRFSS data" / "BRFSS_median_tables_long.csv",
}

THRESHOLD_EVER_CSV_PATHS: dict[str, Path] = {
    "SCCS": ROOT / "SCCS data" / "SCCS_threshold_table_ever_long.csv",
    "MECS": ROOT / "MECS data" / "MECS_threshold_table_ever_long.csv",
    "BRFSS": ROOT / "BRFSS data" / "BRFSS_threshold_table_ever_long.csv",
}

THRESHOLD_TOTAL_CSV_PATHS: dict[str, Path] = {
    "SCCS": ROOT / "SCCS data" / "SCCS_threshold_table_total_long.csv",
    "MECS": ROOT / "MECS data" / "MECS_threshold_table_total_long.csv",
    "BRFSS": ROOT / "BRFSS data" / "BRFSS_threshold_table_total_long.csv",
}

PROBABILITY_CSV_PATHS: dict[str, Path] = {
    "SCCS": ROOT / "SCCS data" / "SCCS_probability_tables_long.csv",
    "MECS": ROOT / "MECS data" / "MECS_probability_tables_long.csv",
    "BRFSS": ROOT / "BRFSS data" / "BRFSS_probability_tables_long.csv",
}

ELIGIBILITY_HEATMAP_HTML_PATHS: dict[str, Path] = {
    "SCCS": ROOT / "SCCS data" / "comp_hm_SCCS_40.html",
    "MECS": ROOT / "MECS data" / "comp_hm_MECS_40.html",
    "BRFSS": ROOT / "BRFSS data" / "comp_hm_BRFSS_40.html",
}

OUTPUTS = {
    "bw": {
        "svg": OUTPUT_DIR / "aggregate_median_table_bw.svg",
        "png": OUTPUT_DIR / "aggregate_median_table_bw.png",
    },
    "full": {
        "svg": OUTPUT_DIR / "aggregate_median_table_full.svg",
        "png": OUTPUT_DIR / "aggregate_median_table_full.png",
    },
    "threshold_ever_bw": {
        "svg": OUTPUT_DIR / "aggregate_threshold_table_ever_bw.svg",
        "png": OUTPUT_DIR / "aggregate_threshold_table_ever_bw.png",
    },
    "threshold_ever_full": {
        "svg": OUTPUT_DIR / "aggregate_threshold_table_ever_full.svg",
        "png": OUTPUT_DIR / "aggregate_threshold_table_ever_full.png",
    },
    "threshold_total_bw": {
        "svg": OUTPUT_DIR / "aggregate_threshold_table_total_bw.svg",
        "png": OUTPUT_DIR / "aggregate_threshold_table_total_bw.png",
    },
    "threshold_total_full": {
        "svg": OUTPUT_DIR / "aggregate_threshold_table_total_full.svg",
        "png": OUTPUT_DIR / "aggregate_threshold_table_total_full.png",
    },
    "probability_p_full": {
        "svg": OUTPUT_DIR / "aggregate_probability_p_full.svg",
        "png": OUTPUT_DIR / "aggregate_probability_p_full.png",
    },
    "probability_2pminus1_full": {
        "svg": OUTPUT_DIR / "aggregate_probability_2pminus1_full.svg",
        "png": OUTPUT_DIR / "aggregate_probability_2pminus1_full.png",
    },
    "figure1": {
        "svg": OUTPUT_DIR / "aggregate_figure1.svg",
        "png": OUTPUT_DIR / "aggregate_figure1.png",
    },
}

COHORT_ORDER = ["SCCS", "MECS", "BRFSS"]
POPULATION_MAP = {
    "Black": "Blacks",
    "White": "Whites",
}

FULL_POPULATION_MAP = {
    "Black": "Blacks",
    "White": "Whites",
    "Black Female": "Black Females",
    "White Female": "White Females",
    "Black Male": "Black Males",
    "White Male": "White Males",
}

MEASURE_COLUMN_MAP = {
    "Pack-Years": "PackYears",
    "Duration": "SmokingYears",
    "Cigarettes Per Day": "CigsPerDay",
}

THRESHOLD_COLUMN_MAP = {
    "≥ 20-Pack Years": "PY20",
    "≥ 30-Pack Years": "PY30",
    "≥ 20-Years": "SY20",
    "≥ 30-Years": "SY30",
}

FIGURE1_MEASURE_MAP = {
    "2021 USPSTF": "Eligible_2021",
    "NCCN (20-Year Duration)": "Eligible_SY20",
}

THRESHOLD_TABLE_NOTE = (
    "Difference calculated as simulated minus observed in percentage points. "
    '"Diff in diff" denotes the difference in these differences between Black and White participants.'
)

MEDIAN_TABLE_NOTE = (
    "Values are reported as median (IQR). Percent difference calculated from simulated and observed "
    "medians as (simulated - observed) / observed x 100."
)

PROBABILITY_HIGHER_NOTE = (
    "Values summarize how often simulated SHG values are higher than observed cohort values.\n"
    "Difference denotes Black minus White within each cohort."
)

PROBABILITY_DIRECTION_NOTE = (
    "Values summarize whether simulated SHG values tend to be higher or lower than observed cohort values, with 0 indicating no\n"
    "directional difference. Positive values indicate higher simulated values and negative values indicate lower simulated values.\n"
    "Difference denotes Black minus White within each cohort."
)

BLACK_RACE_COLOR = "#4292c6"
WHITE_RACE_COLOR = "#084594"
OBSERVED_ORANGE = "orange"
SIMULATED_BLUE = "blue"

# ---------------------------------------------------------------------------
# Loaders and value-formatting helpers
# ---------------------------------------------------------------------------


def load_median_tables() -> dict[str, pd.DataFrame]:
    """Load the finished cohort-specific median table CSVs."""
    loaded: dict[str, pd.DataFrame] = {}

    for cohort, path in MEDIAN_CSV_PATHS.items():
        if not path.exists():
            raise FileNotFoundError(f"Missing expected median table: {path}")
        loaded[cohort] = pd.read_csv(path)

    return loaded


def load_threshold_tables(paths: dict[str, Path]) -> dict[str, pd.DataFrame]:
    """Load the finished cohort-specific threshold table CSVs."""
    loaded: dict[str, pd.DataFrame] = {}

    for cohort, path in paths.items():
        if not path.exists():
            raise FileNotFoundError(f"Missing expected threshold table: {path}")
        loaded[cohort] = pd.read_csv(path)

    return loaded


def load_probability_tables() -> dict[str, pd.DataFrame]:
    """Load the finished cohort-specific probability table CSVs."""
    loaded: dict[str, pd.DataFrame] = {}

    for cohort, path in PROBABILITY_CSV_PATHS.items():
        if not path.exists():
            raise FileNotFoundError(f"Missing expected probability table: {path}")
        loaded[cohort] = pd.read_csv(path)

    return loaded


def load_eligibility_heatmaps() -> dict[str, pd.DataFrame]:
    """Load the eligibility comparison heatmap HTML tables used by Figure 1."""
    loaded: dict[str, pd.DataFrame] = {}

    for cohort, path in ELIGIBILITY_HEATMAP_HTML_PATHS.items():
        if not path.exists():
            raise FileNotFoundError(f"Missing expected eligibility heatmap HTML: {path}")
        loaded[cohort] = parse_simple_html_table(path.read_text(encoding="utf-8", errors="ignore"))

    return loaded


def format_number(value: float | int | None, decimals: int = 1) -> str:
    """Format numeric output with at least one decimal place."""
    if pd.isna(value):
        return ""
    return f"{value:.{decimals}f}"


def format_median_iqr(row: pd.Series, stub: str) -> str:
    """Format one observed or simulated cell as Median (IQR-low-IQR-high)."""
    median = format_number(row[f"{stub}_Median"], decimals=1)
    iqr_low = format_number(round(row[f"{stub}_IQR_Low"]), decimals=0)
    iqr_high = format_number(round(row[f"{stub}_IQR_High"]), decimals=0)
    return f"{median} ({iqr_low}-{iqr_high})"


def format_difference(row: pd.Series, stub: str) -> str:
    """Format one percent-difference cell using the median-only value."""
    return f"{format_number(row[f'{stub}_Median'], decimals=1)}%"


def get_row(df: pd.DataFrame, population: str, source: str) -> pd.Series:
    """Select one row from one cohort median table."""
    match = df.loc[(df["Population"] == population) & (df["Source"] == source)]
    if match.empty:
        raise ValueError(f"Missing row for Population={population}, Source={source}")
    return match.iloc[0]


def format_percent(value: float | int | None, decimals: int = 1) -> str:
    """Format a percentage value with a percent sign."""
    return f"{format_number(value, decimals=decimals)}%"


def format_threshold_observed_or_difference(row: pd.Series, stub: str) -> str:
    """Format observed/difference threshold cells as percentages only."""
    return format_percent(row[f"{stub}_Value"], decimals=1)


def format_threshold_simulated(row: pd.Series, stub: str) -> str:
    """Format simulated threshold cells as value% (low-high)."""
    value = format_percent(row[f"{stub}_Value"], decimals=1)
    low = format_number(round(row[f"{stub}_Lower_CI"]), decimals=0)
    high = format_number(round(row[f"{stub}_Upper_CI"]), decimals=0)
    return f"{value} ({low}-{high})"


def get_probability_row(df: pd.DataFrame, population: str) -> pd.Series:
    """Select one row from one cohort probability table."""
    match = df.loc[df["Population"] == population]
    if match.empty:
        raise ValueError(f"Missing probability row for Population={population}")
    return match.iloc[0]


def build_stat_row(
    loaded: dict[str, pd.DataFrame],
    population_label: str,
    source: str,
    statistic_label: str,
    display_group: str | None = None,
) -> dict[tuple[str, str], str]:
    """Build one content row across all cohorts for one population/source."""
    first_col = statistic_label
    if display_group is not None:
        first_col = f"<strong>{display_group}</strong><br>{statistic_label}"

    row: dict[tuple[str, str], str] = {
        ("", ""): first_col,
    }

    for measure_label, stub in MEASURE_COLUMN_MAP.items():
        for cohort in COHORT_ORDER:
            cohort_row = get_row(loaded[cohort], population_label, source)
            col_name = (measure_label, cohort)

            if source in {"Observed", "Simulated"}:
                row[col_name] = format_median_iqr(cohort_row, stub)
            else:
                row[col_name] = format_difference(cohort_row, stub)

    return row


# ---------------------------------------------------------------------------
# Table builders
# ---------------------------------------------------------------------------

def build_median_table(
    loaded: dict[str, pd.DataFrame],
    population_map: dict[str, str],
) -> list[dict[str, object]]:
    """Build one aggregate median table for the requested subpopulations."""
    rows: list[dict[str, object]] = []

    first_group = next(iter(population_map))

    for display_group, source_population in population_map.items():
        rows.append(
            {
                "row_class": "group-row" if display_group == first_group else "group-row white-group",
                "label_html": f"<strong>{escape(display_group)}</strong>",
                "cells": {
                    (measure, cohort): ""
                    for measure in MEASURE_COLUMN_MAP
                    for cohort in COHORT_ORDER
                },
            }
        )
        rows.append(
            {
                "row_class": "",
                "label_html": "Observed (Median [IQR])",
                "cells": {
                    (measure, cohort): format_median_iqr(
                        get_row(loaded[cohort], source_population, "Observed"),
                        stub,
                    )
                    for measure, stub in MEASURE_COLUMN_MAP.items()
                    for cohort in COHORT_ORDER
                },
            }
        )
        rows.append(
            {
                "row_class": "",
                "label_html": "Simulated (Median [IQR])",
                "cells": {
                    (measure, cohort): format_median_iqr(
                        get_row(loaded[cohort], source_population, "Simulated"),
                        stub,
                    )
                    for measure, stub in MEASURE_COLUMN_MAP.items()
                    for cohort in COHORT_ORDER
                },
            }
        )
        rows.append(
            {
                "row_class": "difference-row",
                "separator_after": display_group != list(population_map.keys())[-1],
                "label_html": "% Difference",
                "cells": {
                    (measure, cohort): format_difference(
                        get_row(loaded[cohort], source_population, "% Difference"),
                        stub,
                    )
                    for measure, stub in MEASURE_COLUMN_MAP.items()
                    for cohort in COHORT_ORDER
                },
            }
        )

    return rows


def build_threshold_table(
    loaded: dict[str, pd.DataFrame],
    population_map: dict[str, str],
) -> list[dict[str, object]]:
    """Build one aggregate threshold table for the requested subpopulations."""
    rows: list[dict[str, object]] = []
    first_group = next(iter(population_map))

    for display_group, source_population in population_map.items():
        rows.append(
            {
                "row_class": "group-row" if display_group == first_group else "group-row white-group",
                "label_html": f"<strong>{escape(display_group)}</strong>",
                "cells": {
                    (measure, cohort): ""
                    for measure in THRESHOLD_COLUMN_MAP
                    for cohort in COHORT_ORDER
                },
            }
        )
        rows.append(
            {
                "row_class": "",
                "label_html": "Observed",
                "cells": {
                    (measure, cohort): format_threshold_observed_or_difference(
                        get_row(loaded[cohort], source_population, "Observed"),
                        stub,
                    )
                    for measure, stub in THRESHOLD_COLUMN_MAP.items()
                    for cohort in COHORT_ORDER
                },
            }
        )
        rows.append(
            {
                "row_class": "",
                "label_html": "Simulated",
                "cells": {
                    (measure, cohort): format_threshold_simulated(
                        get_row(loaded[cohort], source_population, "Simulated"),
                        stub,
                    )
                    for measure, stub in THRESHOLD_COLUMN_MAP.items()
                    for cohort in COHORT_ORDER
                },
            }
        )
        rows.append(
            {
                "row_class": "difference-row",
                "separator_after": True,
                "label_html": "% Difference",
                "cells": {
                    (measure, cohort): format_threshold_observed_or_difference(
                        get_row(loaded[cohort], source_population, "Difference"),
                        stub,
                    )
                    for measure, stub in THRESHOLD_COLUMN_MAP.items()
                    for cohort in COHORT_ORDER
                },
            }
        )

    if list(population_map.keys()) == ["Black", "White"]:
        rows.append(
            {
                "row_class": "diff-in-diff-row",
                "label_html": "Diff in Diff",
                "cells": {
                    (measure, cohort): format_percent(
                        get_row(loaded[cohort], population_map["Black"], "Difference")[f"{stub}_Value"]
                        - get_row(loaded[cohort], population_map["White"], "Difference")[f"{stub}_Value"],
                        decimals=1,
                    )
                    for measure, stub in THRESHOLD_COLUMN_MAP.items()
                    for cohort in COHORT_ORDER
                },
            }
        )

    return rows


def build_threshold_table_full_pairs(loaded: dict[str, pd.DataFrame]) -> list[dict[str, object]]:
    """Build the full aggregate threshold table in paired subgroup blocks."""
    pair_blocks = [
        ("Black", "Blacks", "White", "Whites"),
        ("Black Female", "Black Females", "White Female", "White Females"),
        ("Black Male", "Black Males", "White Male", "White Males"),
    ]

    rows: list[dict[str, object]] = []

    for block_idx, (left_label, left_pop, right_label, right_pop) in enumerate(pair_blocks):
        is_first_group = block_idx == 0

        for display_group, source_population in [(left_label, left_pop), (right_label, right_pop)]:
            rows.append(
                {
                    "row_class": "group-row" if is_first_group and display_group == left_label else "group-row white-group",
                    "label_html": f"<strong>{escape(display_group)}</strong>",
                    "cells": {
                        (measure, cohort): ""
                        for measure in THRESHOLD_COLUMN_MAP
                        for cohort in COHORT_ORDER
                    },
                }
            )
            rows.append(
                {
                    "row_class": "",
                    "label_html": "Observed",
                    "cells": {
                        (measure, cohort): format_threshold_observed_or_difference(
                            get_row(loaded[cohort], source_population, "Observed"),
                            stub,
                        )
                        for measure, stub in THRESHOLD_COLUMN_MAP.items()
                        for cohort in COHORT_ORDER
                    },
                }
            )
            rows.append(
                {
                    "row_class": "",
                    "label_html": "Simulated",
                    "cells": {
                        (measure, cohort): format_threshold_simulated(
                            get_row(loaded[cohort], source_population, "Simulated"),
                            stub,
                        )
                        for measure, stub in THRESHOLD_COLUMN_MAP.items()
                        for cohort in COHORT_ORDER
                    },
                }
            )
            rows.append(
                {
                    "row_class": "difference-row",
                    "separator_after": display_group == right_label,
                    "label_html": "% Difference",
                    "cells": {
                        (measure, cohort): format_threshold_observed_or_difference(
                            get_row(loaded[cohort], source_population, "Difference"),
                            stub,
                        )
                        for measure, stub in THRESHOLD_COLUMN_MAP.items()
                        for cohort in COHORT_ORDER
                    },
                }
            )

        rows.append(
            {
                "row_class": "diff-in-diff-row",
                "separator_after": block_idx < len(pair_blocks) - 1,
                "label_html": "Diff in Diff",
                "cells": {
                    (measure, cohort): format_percent(
                        get_row(loaded[cohort], left_pop, "Difference")[f"{stub}_Value"]
                        - get_row(loaded[cohort], right_pop, "Difference")[f"{stub}_Value"],
                        decimals=1,
                    )
                    for measure, stub in THRESHOLD_COLUMN_MAP.items()
                    for cohort in COHORT_ORDER
                },
            }
        )

    return rows


def build_probability_table_full(
    loaded: dict[str, pd.DataFrame],
    population_map: dict[str, str],
    suffix: str,
) -> list[dict[str, object]]:
    """Build one full aggregate probability table."""
    probability_column_map = {
        "Pack-Years": f"PackYears_{suffix}",
        "Duration": f"SmokingYears_{suffix}",
        "Cigarettes Per Day": f"CigsPerDay_{suffix}",
    }

    rows: list[dict[str, object]] = []

    pair_blocks = [
        ("Black", "Blacks", "White", "Whites"),
        ("Black Female", "Black Females", "White Female", "White Females"),
        ("Black Male", "Black Males", "White Male", "White Males"),
    ]

    for block_idx, (left_label, left_pop, right_label, right_pop) in enumerate(pair_blocks):
        rows.append(
            {
                "row_class": "",
                "label_html": escape(left_label),
                "cells": {
                    (measure, cohort): format_number(
                        get_probability_row(loaded[cohort], left_pop)[column_name],
                        decimals=3,
                    )
                    for measure, column_name in probability_column_map.items()
                    for cohort in COHORT_ORDER
                },
            }
        )
        rows.append(
            {
                "row_class": "",
                "label_html": escape(right_label),
                "cells": {
                    (measure, cohort): format_number(
                        get_probability_row(loaded[cohort], right_pop)[column_name],
                        decimals=3,
                    )
                    for measure, column_name in probability_column_map.items()
                    for cohort in COHORT_ORDER
                },
            }
        )
        rows.append(
            {
                "row_class": "difference-row",
                "separator_after": block_idx < len(pair_blocks) - 1,
                "label_html": "Difference",
                "cells": {
                    (measure, cohort): format_number(
                        get_probability_row(loaded[cohort], left_pop)[column_name]
                        - get_probability_row(loaded[cohort], right_pop)[column_name],
                        decimals=3,
                    )
                    for measure, column_name in probability_column_map.items()
                    for cohort in COHORT_ORDER
                },
            }
        )

    return rows


def parse_simple_html_table(html_text: str) -> pd.DataFrame:
    """Parse the small HTML tables produced by the cohort R scripts."""

    class TableParser(HTMLParser):
        def __init__(self) -> None:
            super().__init__()
            self.in_tr = False
            self.cell_tag: str | None = None
            self.current_cell: list[str] = []
            self.current_row: list[str] = []
            self.headers: list[str] = []
            self.rows: list[list[str]] = []

        def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
            if tag == "tr":
                self.in_tr = True
                self.current_row = []
            elif self.in_tr and tag in {"th", "td"}:
                self.cell_tag = tag
                self.current_cell = []

        def handle_endtag(self, tag: str) -> None:
            if tag in {"th", "td"} and self.cell_tag == tag:
                text = clean_html_cell_text("".join(self.current_cell))
                self.current_row.append(text)
                self.cell_tag = None
                self.current_cell = []
            elif tag == "tr" and self.in_tr:
                if self.current_row:
                    if not self.headers:
                        self.headers = self.current_row
                    else:
                        self.rows.append(self.current_row)
                self.in_tr = False
                self.current_row = []

        def handle_data(self, data: str) -> None:
            if self.cell_tag is not None:
                self.current_cell.append(data)

    parser = TableParser()
    parser.feed(html_text)

    if not parser.headers:
        raise ValueError("No header row found in HTML table.")

    valid_rows = [
        row for row in parser.rows
        if len(row) == len(parser.headers) and any(cell.strip() for cell in row)
    ]

    return pd.DataFrame(valid_rows, columns=parser.headers)


def clean_html_cell_text(text: str) -> str:
    """Strip inline HTML tags and normalize whitespace."""
    cleaned = re.sub(r"<[^>]+>", "", text)
    cleaned = cleaned.replace("&nbsp;", " ")
    return " ".join(cleaned.split())


def parse_signed_percent(text: str) -> float:
    """Convert a signed percentage string like '+31.1%' to float."""
    cleaned = text.strip().replace("%", "")
    return float(cleaned)


def build_figure1_data(loaded: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build the top-panel bar data and bottom table data for Figure 1."""
    plot_records: list[dict[str, object]] = []
    table_records: list[dict[str, object]] = []

    for measure_label, source_measure in FIGURE1_MEASURE_MAP.items():
        for cohort in COHORT_ORDER:
            df = loaded[cohort]
            measure_col = df.columns[0]
            row = df.loc[df[measure_col] == source_measure]
            if row.empty:
                raise ValueError(f"Missing {source_measure} row in {cohort} eligibility table.")
            row = row.iloc[0]

            black_val = parse_signed_percent(str(row["Black"]))
            white_val = parse_signed_percent(str(row["White"]))

            plot_records.extend(
                [
                    {
                        "Measure": measure_label,
                        "Cohort": cohort,
                        "Race": "Black",
                        "Value": black_val,
                    },
                    {
                        "Measure": measure_label,
                        "Cohort": cohort,
                        "Race": "White",
                        "Value": white_val,
                    },
                ]
            )

            table_records.append(
                {
                    "Population": "Black",
                    "Measure": measure_label,
                    "Cohort": cohort,
                    "Value": black_val,
                }
            )
            table_records.append(
                {
                    "Population": "White",
                    "Measure": measure_label,
                    "Cohort": cohort,
                    "Value": white_val,
                }
            )
            table_records.append(
                {
                    "Population": "Difference",
                    "Measure": measure_label,
                    "Cohort": cohort,
                    "Value": black_val - white_val,
                }
            )

    return pd.DataFrame(plot_records), pd.DataFrame(table_records)


def render_html_table(rows: list[dict[str, object]]) -> str:
    """Render a manual HTML table.

    This HTML path is retained as a readable prototype/output reference.
    The main paper-ready deliverables are the SVG and PNG exports.
    """
    header_row_1 = (
        "<tr>"
        "<th></th>"
        f"<th colspan='3'>Pack-Years</th>"
        f"<th colspan='3'>Duration</th>"
        f"<th colspan='3'>Cigarettes Per Day</th>"
        "</tr>"
    )

    header_row_2 = (
        "<tr>"
        "<th></th>"
        + "".join(f"<th>{escape(cohort)}</th>" for _measure in THRESHOLD_COLUMN_MAP if False)
        + "</tr>"
    )
    # The HTML export is no longer a primary deliverable, but keep it structurally aligned.
    measure_keys = list(rows[0]["cells"].keys()) if rows else []
    header_row_2 = (
        "<tr>"
        "<th></th>"
        + "".join(f"<th>{escape(cohort)}</th>" for _measure in sorted({m for m, _c in measure_keys}, key=lambda x: list(dict.fromkeys([m for m, _ in measure_keys])).index(x)) for cohort in COHORT_ORDER)
        + "</tr>"
    )

    body_rows: list[str] = []
    for row in rows:
        label_html = row["label_html"]
        cells = row["cells"]
        row_class = row.get("row_class", "")
        class_attr = f" class='{row_class}'" if row_class else ""
        row_html = f"<tr{class_attr}>"
        row_html += f"<td class='label'>{label_html}</td>"
        for measure in MEASURE_COLUMN_MAP:
            for cohort in COHORT_ORDER:
                row_html += f"<td>{escape(str(cells[(measure, cohort)]))}</td>"
        row_html += "</tr>"
        body_rows.append(row_html)

    return (
        "<table>"
        "<thead>"
        f"{header_row_1}"
        f"{header_row_2}"
        "</thead>"
        "<tbody>"
        + "".join(body_rows)
        + "</tbody>"
        "</table>"
    )


def save_styled_html(rows: list[dict[str, object]], output_path: Path) -> None:
    """Save a simple manual HTML version of the table for quick inspection."""
    table_html = render_html_table(rows)

    full_html = f"""
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Aggregate Median Table</title>
  <style>
    body {{
      font-family: Arial, sans-serif;
      margin: 24px;
    }}
    table {{
      border-collapse: collapse;
      font-size: 13px;
    }}
    th, td {{
      padding: 6px;
      text-align: center;
      vertical-align: middle;
      background: white;
      border: none;
    }}
    thead th {{
      font-weight: bold;
      background: white;
    }}
    thead tr:first-child th {{
      font-size: 15px;
    }}
    thead tr:nth-child(2) th {{
      font-size: 13px;
      font-weight: normal;
    }}
    td.label {{
      text-align: left;
      white-space: nowrap;
    }}
    thead tr:first-child th {{
      border-bottom: 1px solid black;
    }}
    thead tr:first-child th:nth-child(1) {{
      border-right: 1px solid black;
    }}
    thead tr:nth-child(2) th:nth-child(1),
    tbody td:nth-child(1) {{
      border-right: 1px solid black;
    }}
    thead tr:nth-child(2) th:nth-child(4),
    tbody td:nth-child(4) {{
      border-right: 1px solid black;
    }}
    thead tr:nth-child(2) th:nth-child(7),
    tbody td:nth-child(7) {{
      border-right: 1px solid black;
    }}
    thead tr:first-child th:nth-child(2),
    thead tr:first-child th:nth-child(3) {{
      border-right: 1px solid black;
    }}
    tbody tr.difference-row:first-of-type td {{
      border-bottom: 1px solid black;
    }}
    tbody tr.white-group td {{
      border-top: 1px solid black;
    }}
  </style>
</head>
<body>
{table_html}
</body>
</html>
"""

    output_path.write_text(full_html, encoding="utf-8")


def make_border_xml(side: str) -> str:
    return f'<w:{side} w:val="single" w:sz="8" w:space="0" w:color="000000"/>'


def make_cell_xml(
    text: str = "",
    *,
    width: int,
    bold: bool = False,
    align: str = "center",
    colspan: int = 1,
    right: bool = False,
    bottom: bool = False,
) -> str:
    """Create one Word table cell."""
    tc_pr = [f'<w:tcW w:w="{width}" w:type="dxa"/>']

    if colspan > 1:
        tc_pr.append(f'<w:gridSpan w:val="{colspan}"/>')

    borders: list[str] = []
    if right:
        borders.append(make_border_xml("right"))
    if bottom:
        borders.append(make_border_xml("bottom"))
    if borders:
        tc_pr.append(f"<w:tcBorders>{''.join(borders)}</w:tcBorders>")

    jc = "left" if align == "left" else "center"
    runs: list[str] = []
    for i, part in enumerate(text.split("<br>")):
        run_props = "<w:rPr><w:b/></w:rPr>" if bold else ""
        runs.append(f"<w:r>{run_props}<w:t xml:space=\"preserve\">{xml_escape(part)}</w:t></w:r>")
        if i < len(text.split("<br>")) - 1:
            runs.append("<w:r><w:br/></w:r>")

    paragraph = (
        f"<w:p><w:pPr><w:jc w:val=\"{jc}\"/></w:pPr>"
        f"{''.join(runs)}"
        f"</w:p>"
    )

    return f"<w:tc><w:tcPr>{''.join(tc_pr)}</w:tcPr>{paragraph}</w:tc>"


def build_docx_table_xml(rows: list[dict[str, object]]) -> str:
    """Build a minimal WordprocessingML table.

    This DOCX path is kept as a fallback export route, but the current
    manuscript-facing outputs are SVG and PNG.
    """
    label_w = 3600
    data_w = 1200
    grid_cols = [label_w] + [data_w] * 9

    parts: list[str] = [
        "<w:tbl>",
        (
            "<w:tblPr>"
            "<w:tblLayout w:type=\"fixed\"/>"
            "<w:tblW w:w=\"0\" w:type=\"auto\"/>"
            "</w:tblPr>"
        ),
        "<w:tblGrid>" + "".join(f'<w:gridCol w:w="{w}"/>' for w in grid_cols) + "</w:tblGrid>",
    ]

    # Header row 1
    parts.append(
        "<w:tr>"
        + make_cell_xml("", width=label_w, right=True, bottom=True)
        + make_cell_xml("Pack-Years", width=data_w * 3, colspan=3, right=True, bottom=True)
        + make_cell_xml("Duration", width=data_w * 3, colspan=3, right=True, bottom=True)
        + make_cell_xml("Cigarettes Per Day", width=data_w * 3, colspan=3, bottom=True)
        + "</w:tr>"
    )

    # Header row 2
    header_cells = [make_cell_xml("", width=label_w, right=True)]
    for idx, cohort in enumerate(COHORT_ORDER, start=1):
        header_cells.append(make_cell_xml(cohort, width=data_w))
    for idx, cohort in enumerate(COHORT_ORDER, start=1):
        header_cells.append(make_cell_xml(cohort, width=data_w, right=(idx == 3)))
    for idx, cohort in enumerate(COHORT_ORDER, start=1):
        header_cells.append(make_cell_xml(cohort, width=data_w, right=(idx == 3)))
    # fix right borders after col 4 and 7 overall
    header_cells = [
        make_cell_xml("", width=label_w, right=True),
        make_cell_xml("SCCS", width=data_w),
        make_cell_xml("MECS", width=data_w),
        make_cell_xml("BRFSS", width=data_w, right=True),
        make_cell_xml("SCCS", width=data_w),
        make_cell_xml("MECS", width=data_w),
        make_cell_xml("BRFSS", width=data_w, right=True),
        make_cell_xml("SCCS", width=data_w),
        make_cell_xml("MECS", width=data_w),
        make_cell_xml("BRFSS", width=data_w),
    ]
    parts.append("<w:tr>" + "".join(header_cells) + "</w:tr>")

    # Body rows
    first_difference_idx = next(
        (idx for idx, row in enumerate(rows) if row.get("row_class", "") == "difference-row"),
        None,
    )

    for idx, row in enumerate(rows):
        label_text = str(row["label_html"]).replace("<strong>", "").replace("</strong>", "")
        cells = row["cells"]
        is_group_row = row.get("row_class", "") in {"group-row", "group-row white-group"}
        add_bottom = idx == first_difference_idx

        body_cells = [
            make_cell_xml(
                label_text,
                width=label_w,
                bold=is_group_row,
                align="left",
                right=True,
                bottom=add_bottom,
            )
        ]

        data_order = [(measure, cohort) for measure in MEASURE_COLUMN_MAP for cohort in COHORT_ORDER]
        for cell_idx, key in enumerate(data_order, start=1):
            body_cells.append(
                make_cell_xml(
                    str(cells[key]),
                    width=data_w,
                    right=(cell_idx in {3, 6}),
                    bottom=add_bottom,
                )
            )

        parts.append("<w:tr>" + "".join(body_cells) + "</w:tr>")

    parts.append("</w:tbl>")
    return "".join(parts)


def save_docx(rows: list[dict[str, object]], output_path: Path) -> None:
    """Write a minimal .docx file containing the aggregate median table."""
    created = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    table_xml = build_docx_table_xml(rows)

    document_xml = f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:document xmlns:wpc="http://schemas.microsoft.com/office/word/2010/wordprocessingCanvas"
 xmlns:mc="http://schemas.openxmlformats.org/markup-compatibility/2006"
 xmlns:o="urn:schemas-microsoft-com:office:office"
 xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"
 xmlns:m="http://schemas.openxmlformats.org/officeDocument/2006/math"
 xmlns:v="urn:schemas-microsoft-com:vml"
 xmlns:wp14="http://schemas.microsoft.com/office/word/2010/wordprocessingDrawing"
 xmlns:wp="http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing"
 xmlns:w10="urn:schemas-microsoft-com:office:word"
 xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"
 xmlns:w14="http://schemas.microsoft.com/office/word/2010/wordml"
 xmlns:wpg="http://schemas.microsoft.com/office/word/2010/wordprocessingGroup"
 xmlns:wpi="http://schemas.microsoft.com/office/word/2010/wordprocessingInk"
 xmlns:wne="http://schemas.microsoft.com/office/word/2006/wordml"
 xmlns:wps="http://schemas.microsoft.com/office/word/2010/wordprocessingShape"
 mc:Ignorable="w14 wp14">
  <w:body>
    {table_xml}
    <w:sectPr>
      <w:pgSz w:w="12240" w:h="15840"/>
      <w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" w:header="708" w:footer="708" w:gutter="0"/>
    </w:sectPr>
  </w:body>
</w:document>
"""

    content_types_xml = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>
  <Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>
  <Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>
</Types>
"""

    rels_xml = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>
  <Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>
</Relationships>
"""

    document_rels_xml = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"/>
"""

    core_xml = f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties"
 xmlns:dc="http://purl.org/dc/elements/1.1/"
 xmlns:dcterms="http://purl.org/dc/terms/"
 xmlns:dcmitype="http://purl.org/dc/dcmitype/"
 xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <dc:title>Aggregate Median Table</dc:title>
  <dc:creator>Hesam Mahmoudi</dc:creator>
  <cp:lastModifiedBy>Hesam Mahmoudi</cp:lastModifiedBy>
  <dcterms:created xsi:type="dcterms:W3CDTF">{created}</dcterms:created>
  <dcterms:modified xsi:type="dcterms:W3CDTF">{created}</dcterms:modified>
</cp:coreProperties>
"""

    app_xml = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties"
 xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">
  <Application>Microsoft Office Word</Application>
</Properties>
"""

    with zipfile.ZipFile(output_path, "w", compression=zipfile.ZIP_DEFLATED) as docx:
        docx.writestr("[Content_Types].xml", content_types_xml)
        docx.writestr("_rels/.rels", rels_xml)
        docx.writestr("word/document.xml", document_xml)
        docx.writestr("word/_rels/document.xml.rels", document_rels_xml)
        docx.writestr("docProps/core.xml", core_xml)
        docx.writestr("docProps/app.xml", app_xml)


def strip_html_text(text: str) -> list[str]:
    """Remove the limited inline markup used in labels and return linewise text."""
    cleaned = (
        text.replace("<strong>", "")
        .replace("</strong>", "")
        .replace("<br>", "\n")
    )
    return cleaned.split("\n")


def compute_column_widths(rows: list[dict[str, object]]) -> list[float]:
    """Compute content-aware widths for the label column and 9 data columns."""
    label_lines = [""]
    for row in rows:
        label_lines.extend(strip_html_text(str(row["label_html"])))

    label_chars = max(len(line) for line in label_lines if line is not None)
    label_width = min(max(0.085 * label_chars + 0.35, 1.9), 3.0)

    data_order = [(measure, cohort) for measure in MEASURE_COLUMN_MAP for cohort in COHORT_ORDER]
    data_widths: list[float] = []

    for measure, cohort in data_order:
        values = [cohort]
        values.extend(str(row["cells"][(measure, cohort)]) for row in rows)
        max_chars = max(len(value) for value in values)
        width = min(max(0.065 * max_chars + 0.30, 0.82), 1.35)
        data_widths.append(width)

    # Ensure each 3-column measure block is wide enough for its group title.
    for block_idx, measure in enumerate(MEASURE_COLUMN_MAP):
        start = block_idx * 3
        end = start + 3
        current_sum = sum(data_widths[start:end])
        min_sum = min(max(0.05 * len(measure) + 0.95, 2.65), 3.9)
        if current_sum < min_sum:
            extra = (min_sum - current_sum) / 3
            for idx in range(start, end):
                data_widths[idx] += extra

    return [label_width] + data_widths


def compute_column_widths_generic(
    rows: list[dict[str, object]],
    measure_map: dict[str, str],
) -> list[float]:
    """Compute content-aware widths for a generic grouped table."""
    label_lines = [""]
    for row in rows:
        label_lines.extend(strip_html_text(str(row["label_html"])))

    label_chars = max(len(line) for line in label_lines if line is not None)
    if measure_map == THRESHOLD_COLUMN_MAP:
        label_width = min(max(0.062 * label_chars + 0.12, 1.30), 2.00)
    else:
        label_width = min(max(0.085 * label_chars + 0.35, 1.9), 3.0)

    data_order = [(measure, cohort) for measure in measure_map for cohort in COHORT_ORDER]
    data_widths: list[float] = []

    for measure, cohort in data_order:
        values = [cohort]
        values.extend(str(row["cells"][(measure, cohort)]) for row in rows)
        max_chars = max(len(value) for value in values)
        width = min(max(0.065 * max_chars + 0.30, 0.82), 1.35)
        data_widths.append(width)

    for block_idx, measure in enumerate(measure_map):
        start = block_idx * 3
        end = start + 3
        current_sum = sum(data_widths[start:end])
        min_sum = min(max(0.05 * len(measure) + 0.95, 2.65), 3.9)
        if current_sum < min_sum:
            extra = (min_sum - current_sum) / 3
            for idx in range(start, end):
                data_widths[idx] += extra

    return [label_width] + data_widths


# ---------------------------------------------------------------------------
# Figure rendering
# ---------------------------------------------------------------------------

def save_table_figure(rows: list[dict[str, object]], svg_path: Path, png_path: Path) -> None:
    """Render a grouped aggregate table using the default 3-measure layout."""
    save_table_figure_generic(
        rows,
        MEASURE_COLUMN_MAP,
        svg_path,
        png_path,
        note_text=MEDIAN_TABLE_NOTE,
        note_wrap_width=200,
    )


def save_table_figure_generic(
    rows: list[dict[str, object]],
    measure_map: dict[str, str],
    svg_path: Path,
    png_path: Path,
    note_text: str | None = None,
    note_wrap_width: int | None = None,
    label_width_override: float | None = None,
) -> None:
    """Render a generic grouped aggregate table as SVG and PNG.

    The figure layout is shared across median, threshold, and probability
    tables. Widths are computed from the actual cell contents so the tables can
    be reused across outputs without hand-tuning each column.
    """
    plt.rcParams["font.family"] = "Arial"
    col_widths = compute_column_widths_generic(rows, measure_map)
    if label_width_override is not None:
        # Probability tables read better with a narrower label stub than the
        # shared generic minimum that works for median tables.
        col_widths[0] = label_width_override
    row_heights = [0.76, 0.62] + [0.60] * len(rows)
    # Keep note text inside the drawn table width so savefig does not expand or
    # visually compress the table area to accommodate overflow.
    wrap_width = note_wrap_width if note_wrap_width is not None else max(60, int(sum(col_widths) * 11.5))
    note_lines: list[str] = []
    if note_text:
        raw_lines = note_text.splitlines()
        if len(raw_lines) > 1:
            # Respect manual line breaks when the note wording has already been
            # tuned to the existing table width.
            note_lines = raw_lines
        else:
            note_lines = wrap(note_text, width=wrap_width)
    note_line_height = 0.34
    note_top_gap = 0.18 if note_lines else 0.0
    note_bottom_gap = 0.10 if note_lines else 0.0

    total_width = sum(col_widths)
    table_height = sum(row_heights)
    total_height = table_height + note_top_gap + len(note_lines) * note_line_height + note_bottom_gap

    fig, ax = plt.subplots(figsize=(total_width * 0.76, total_height * 0.38), dpi=300)
    ax.set_xlim(0, total_width)
    ax.set_ylim(total_height, 0)
    ax.axis("off")
    ax.set_frame_on(False)
    ax.patch.set_visible(False)

    x_edges = [0]
    for width in col_widths:
        x_edges.append(x_edges[-1] + width)

    y_edges = [0]
    for height in row_heights:
        y_edges.append(y_edges[-1] + height)

    # Header row 1
    for block_idx, measure in enumerate(measure_map):
        start_col = 1 + block_idx * 3
        end_col = start_col + 3
        ax.text(
            (x_edges[start_col] + x_edges[end_col]) / 2,
            (y_edges[0] + y_edges[1]) / 2,
            measure,
            ha="center",
            va="center",
            fontsize=11,
            fontweight="bold",
        )

    # Header row 2
    cohorts = COHORT_ORDER * len(measure_map)
    for idx, cohort in enumerate(cohorts, start=1):
        x_center = (x_edges[idx] + x_edges[idx + 1]) / 2
        y_center = (y_edges[1] + y_edges[2]) / 2
        ax.text(x_center, y_center, cohort, ha="center", va="center", fontsize=10)

    # Body rows
    for i, row in enumerate(rows, start=2):
        y_center = (y_edges[i] + y_edges[i + 1]) / 2
        label_text = str(row["label_html"]).replace("<strong>", "").replace("</strong>", "")
        label_lines = label_text.split("<br>")
        is_group = row.get("row_class", "") in {"group-row", "group-row white-group"}

        if len(label_lines) == 2:
            ax.text(x_edges[0] + 0.08, y_center - 0.11, label_lines[0],
                    ha="left", va="center", fontsize=10.5, fontweight="bold")
            ax.text(x_edges[0] + 0.08, y_center + 0.11, label_lines[1],
                    ha="left", va="center", fontsize=9.2)
        else:
            ax.text(x_edges[0] + 0.08, y_center, label_lines[0],
                    ha="left", va="center", fontsize=9.4,
                    fontweight="bold" if is_group else "normal")

        col_idx = 1
        for measure in measure_map:
            for cohort in COHORT_ORDER:
                value = str(row["cells"][(measure, cohort)])
                x_center = (x_edges[col_idx] + x_edges[col_idx + 1]) / 2
                ax.text(x_center, y_center, value, ha="center", va="center", fontsize=9.2)
                col_idx += 1

    # Vertical dividers
    divider_cols = [1] + [1 + 3 * i for i in range(1, len(measure_map))]
    for divider_col in divider_cols:
        x = x_edges[divider_col]
        ax.plot([x, x], [0, table_height], color="black", linewidth=0.8)

    # Horizontal dividers
    ax.plot([0, total_width], [0, 0], color="black", linewidth=0.8)
    ax.plot([0, total_width], [y_edges[1], y_edges[1]], color="black", linewidth=0.8)
    for idx, row in enumerate(rows):
        if row.get("separator_after", False):
            ax.plot([0, total_width], [y_edges[idx + 3], y_edges[idx + 3]], color="black", linewidth=0.8)

    if note_lines:
        bottom_rule_y = table_height
        ax.plot([0, total_width], [bottom_rule_y, bottom_rule_y], color="black", linewidth=0.8)
        note_y = bottom_rule_y + note_top_gap + note_line_height / 2
        for line in note_lines:
            ax.text(
                0,
                note_y,
                line,
                ha="left",
                va="center",
                fontsize=8.6,
            )
            note_y += note_line_height

    # Fixed margins plus tight cropping keep the table geometry stable while
    # trimming excess canvas around the exported figure.
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    fig.savefig(svg_path, format="svg", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(png_path, format="png", dpi=300, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def save_figure1(plot_df: pd.DataFrame, table_df: pd.DataFrame, svg_path: Path, png_path: Path) -> None:
    """Render the manuscript combined bar-chart/table summary figure.

    The underlying structure follows the supervisor's original layout, while
    the colors and labels are aligned to the rest of the paper:
    orange/blue remain reserved for observed-vs-simulated interpretation, and
    race colors follow the Black/White subgroup palette used in the forest
    plots.
    """
    plt.rcParams["font.family"] = "Arial"

    fig = plt.figure(figsize=(7.65, 5.2), dpi=300)
    gs = fig.add_gridspec(
        nrows=2,
        ncols=2,
        width_ratios=[14.0, 2.8],
        height_ratios=[3.8, 1.15],
        hspace=0.0,
        wspace=0.03,
    )

    ax = fig.add_subplot(gs[0, 0])
    ax_annot = fig.add_subplot(gs[0, 1])
    ax_table = fig.add_subplot(gs[1, 0])
    ax_blank = fig.add_subplot(gs[1, 1])
    ax_blank.axis("off")
    ax_annot.set_zorder(5)
    ax_annot.patch.set_alpha(0)

    measure_order = list(FIGURE1_MEASURE_MAP.keys())
    x_left = 0.0
    x_right = 4.8
    unit = (x_right - x_left) / 6.0
    cohort_centers = [x_left + (i + 0.5) * unit for i in range(6)]
    cohort_positions = {
        measure_order[0]: {
            "SCCS": cohort_centers[0],
            "MECS": cohort_centers[1],
            "BRFSS": cohort_centers[2],
        },
        measure_order[1]: {
            "SCCS": cohort_centers[3],
            "MECS": cohort_centers[4],
            "BRFSS": cohort_centers[5],
        },
    }
    bar_width = 0.225
    for measure in measure_order:
        for cohort in COHORT_ORDER:
            base_x = cohort_positions[measure][cohort]
            for race, offset in [("Black", -bar_width / 2), ("White", bar_width / 2)]:
                value = float(
                    plot_df.loc[
                        (plot_df["Measure"] == measure)
                        & (plot_df["Cohort"] == cohort)
                        & (plot_df["Race"] == race),
                        "Value",
                    ].iloc[0]
                )
                ax.bar(
                    base_x + offset,
                    value,
                    width=bar_width * 0.92,
                    color="white" if race == "Black" else WHITE_RACE_COLOR,
                    edgecolor=BLACK_RACE_COLOR if race == "Black" else WHITE_RACE_COLOR,
                    linewidth=1.0 if race == "Black" else 0.6,
                    hatch="///" if race == "Black" else None,
                    label=race if (measure == measure_order[0] and cohort == "SCCS") else None,
                )
                label_y = value + 2.2 if value >= 0 else value - 3.2
                ax.text(
                    base_x + offset,
                    label_y,
                    f"{int(round(value))}",
                    ha="center",
                    va="bottom" if value >= 0 else "top",
                    fontsize=10,
                    fontweight="bold" if abs(value) >= 80 else "normal",
                )

    ax.axhline(0, color="#f28e2b", linewidth=2.4)
    section_divider_x = x_left + 3 * unit
    plot_x0, plot_x1 = -0.03, x_right
    ax.plot([section_divider_x, section_divider_x], [-36, 116], color="black", linewidth=1.0, linestyle=(0, (4, 3)))
    ax.text(
        plot_x1 + 0.12,
        0,
        "Observed data",
        color="#f28e2b",
        ha="left",
        va="center",
        fontsize=9.6,
        fontweight="bold",
        clip_on=False,
    )
    ax.set_ylabel("Percent Difference", fontsize=11)
    ax.set_xlim(plot_x0, plot_x1)
    ax.set_ylim(-45, 125)
    ax.set_xticks([])
    ax.set_yticks([-40, -20, 0, 20, 40, 60, 80, 100, 120])
    ax.tick_params(axis="y", length=0, labelsize=10)
    for spine in ["top", "right", "left", "bottom"]:
        ax.spines[spine].set_visible(False)

    ax.legend(
        loc="upper right",
        bbox_to_anchor=(0.96, 0.90),
        frameon=False,
        ncol=2,
        fontsize=10.5,
        handlelength=0.9,
        handletextpad=0.3,
        columnspacing=1.1,
    )

    ax_annot.set_xlim(0, 1)
    ax_annot.set_ylim(0, 1)
    ax_annot.axis("off")
    zero_frac = (0 - (-45)) / (125 - (-45))
    upper_start = zero_frac + 0.03
    lower_start = zero_frac - 0.03
    upper_end = 0.80
    lower_end = -0.02
    ax_annot.arrow(
        0.04, upper_start, 0, upper_end - upper_start,
        width=0.014,
        head_width=0.077,
        head_length=0.028,
        color="black",
        length_includes_head=True,
        clip_on=False,
        zorder=5,
    )
    ax_annot.arrow(
        0.04, lower_start, 0, lower_end - lower_start,
        width=0.014,
        head_width=0.077,
        head_length=0.028,
        color="black",
        length_includes_head=True,
        clip_on=False,
        zorder=5,
    )
    upper_center = (upper_start + upper_end) / 2
    lower_center = (lower_start + lower_end) / 2
    line_step = 0.040

    upper_lines = [
        ("SHG", "normal", upper_center + 2 * line_step),
        ("overestimates", "bold", upper_center + 1 * line_step),
        ("eligibility", "normal", upper_center),
        ("compared to", "normal", upper_center - 1 * line_step),
        ("observed data", "normal", upper_center - 2 * line_step),
    ]
    lower_lines = [
        ("SHG", "normal", lower_center + 2 * line_step),
        ("underestimates", "bold", lower_center + 1 * line_step),
        ("eligibility", "normal", lower_center),
        ("compared to", "normal", lower_center - 1 * line_step),
        ("observed data", "normal", lower_center - 2 * line_step),
    ]

    for text, weight, ypos in upper_lines:
        ax_annot.text(
            0.11, ypos, text,
            color="black",
            ha="left",
            va="center",
            fontsize=9.6,
            fontweight=weight,
        )
    for text, weight, ypos in lower_lines:
        ax_annot.text(
            0.11, ypos, text,
            color="black",
            ha="left",
            va="center",
            fontsize=9.6,
            fontweight=weight,
        )

    ax_table.axis("off")
    table_rows = ["Black", "White", "Difference"]
    row_heights = [0.46, 0.46, 0.52, 0.52, 0.52]
    total_height = sum(row_heights)
    ax_table.set_xlim(0, 1)
    ax_table.set_ylim(total_height, 0)
    mixed_transform = transforms.blended_transform_factory(ax_table.transAxes, ax_table.transData)

    y_edges = [0]
    for height in row_heights:
        y_edges.append(y_edges[-1] + height)

    def x_to_frac(x: float) -> float:
        return (x - plot_x0) / (plot_x1 - plot_x0)

    cohort_fracs = [x_to_frac(x) for x in cohort_centers]
    section_frac = x_to_frac(section_divider_x)
    left_boundary_frac = x_to_frac(x_left)
    right_boundary_frac = x_to_frac(x_right)
    label_col_left_frac = -0.155
    marker_x = -0.132
    marker_w = 0.080
    marker_text_x = -0.106
    diff_text_x = marker_x
    cohort_boundaries = [
        x_to_frac(x_left + 1 * unit),
        x_to_frac(x_left + 2 * unit),
        x_to_frac(x_left + 4 * unit),
        x_to_frac(x_left + 5 * unit),
    ]

    for x_center, cohort in zip(cohort_centers, COHORT_ORDER + COHORT_ORDER):
        ax_table.text(
            x_to_frac(x_center),
            (y_edges[0] + y_edges[1]) / 2,
            cohort,
            ha="center",
            va="center",
            fontsize=10.5,
            transform=mixed_transform,
        )

    merged_spans = [
        ("2021 USPSTF", left_boundary_frac, section_frac),
        ("NCCN (20-Year Duration)", section_frac, right_boundary_frac),
    ]
    for measure_name, x0, x1 in merged_spans:
        ax_table.text(
            (x0 + x1) / 2,
            (y_edges[1] + y_edges[2]) / 2,
            measure_name,
            ha="center",
            va="center",
            fontsize=10.5,
            transform=mixed_transform,
        )

    marker_specs = {
        "Black": {"facecolor": "white", "edgecolor": BLACK_RACE_COLOR, "hatch": "///"},
        "White": {"facecolor": WHITE_RACE_COLOR, "edgecolor": WHITE_RACE_COLOR, "hatch": None},
    }

    for row_idx, population in enumerate(table_rows, start=2):
        y_center = (y_edges[row_idx] + y_edges[row_idx + 1]) / 2
        if population in marker_specs:
            spec = marker_specs[population]
            marker_h = 0.18
            rect = plt.Rectangle(
                (marker_x, y_center - marker_h / 2),
                marker_w / (plot_x1 - plot_x0),
                marker_h,
                facecolor=spec["facecolor"],
                edgecolor=spec["edgecolor"],
                linewidth=0.7,
                hatch=spec["hatch"],
                transform=mixed_transform,
                clip_on=False,
            )
            ax_table.add_patch(rect)
            text_x = marker_text_x
        else:
            text_x = diff_text_x

        ax_table.text(
            text_x,
            y_center,
            population,
            ha="left",
            va="center",
            fontsize=10.5,
            fontweight="bold" if population in {"Black", "White"} else "normal",
            transform=mixed_transform,
            clip_on=False,
        )

        col_idx = 0
        for measure in measure_order:
            for cohort in COHORT_ORDER:
                value = float(
                    table_df.loc[
                        (table_df["Population"] == population)
                        & (table_df["Measure"] == measure)
                        & (table_df["Cohort"] == cohort),
                        "Value",
                    ].iloc[0]
                )
                x_center = cohort_fracs[col_idx]
                ax_table.text(
                    x_center,
                    y_center,
                    f"{value:.2f}",
                    ha="center",
                    va="center",
                    fontsize=10.1,
                    transform=mixed_transform,
                )
                col_idx += 1

    line_color = "#d0d0d0"
    # Horizontal rules. Keep the upper-left header area visually open.
    ax_table.plot([left_boundary_frac, right_boundary_frac], [y_edges[0], y_edges[0]], color=line_color, linewidth=0.8, transform=mixed_transform, clip_on=False)
    ax_table.plot([left_boundary_frac, right_boundary_frac], [y_edges[1], y_edges[1]], color=line_color, linewidth=0.8, transform=mixed_transform, clip_on=False)
    for y in y_edges[2:]:
        ax_table.plot([label_col_left_frac, right_boundary_frac], [y, y], color=line_color, linewidth=0.8, transform=mixed_transform, clip_on=False)

    # Outer/major vertical rules go full height.
    for x in [left_boundary_frac, section_frac, right_boundary_frac]:
        ax_table.plot([x, x], [y_edges[0], y_edges[-1]], color=line_color, linewidth=0.8, transform=mixed_transform, clip_on=False)
    ax_table.plot([label_col_left_frac, label_col_left_frac], [y_edges[2], y_edges[-1]], color=line_color, linewidth=0.8, transform=mixed_transform, clip_on=False)

    # Cohort separators skip the merged measure-name row.
    for x in cohort_boundaries:
        ax_table.plot([x, x], [y_edges[0], y_edges[1]], color=line_color, linewidth=0.8, transform=mixed_transform, clip_on=False)
        ax_table.plot([x, x], [y_edges[2], y_edges[-1]], color=line_color, linewidth=0.8, transform=mixed_transform, clip_on=False)

    fig.subplots_adjust(left=0.115, right=0.985, top=0.98, bottom=0.06)
    fig.savefig(svg_path, format="svg", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(png_path, format="png", bbox_inches="tight", pad_inches=0.02, dpi=300)
    plt.close(fig)


def main() -> None:
    """Build all current aggregate table figures from the cohort CSV exports."""
    median_tables = load_median_tables()
    threshold_ever_tables = load_threshold_tables(THRESHOLD_EVER_CSV_PATHS)
    threshold_total_tables = load_threshold_tables(THRESHOLD_TOTAL_CSV_PATHS)
    probability_tables = load_probability_tables()
    eligibility_heatmaps = load_eligibility_heatmaps()
    tables = {
        "bw": build_median_table(median_tables, POPULATION_MAP),
        "full": build_median_table(median_tables, FULL_POPULATION_MAP),
    }

    for key, table_rows in tables.items():
        save_table_figure(table_rows, OUTPUTS[key]["svg"], OUTPUTS[key]["png"])
        print()
        print(f"Saved: {OUTPUTS[key]['svg']}")
        print(f"Saved: {OUTPUTS[key]['png']}")

    threshold_bw = build_threshold_table(threshold_ever_tables, POPULATION_MAP)
    save_table_figure_generic(
        threshold_bw,
        THRESHOLD_COLUMN_MAP,
        OUTPUTS["threshold_ever_bw"]["svg"],
        OUTPUTS["threshold_ever_bw"]["png"],
        note_text=THRESHOLD_TABLE_NOTE,
    )
    print()
    print(f"Saved: {OUTPUTS['threshold_ever_bw']['svg']}")
    print(f"Saved: {OUTPUTS['threshold_ever_bw']['png']}")

    threshold_full = build_threshold_table_full_pairs(threshold_ever_tables)
    save_table_figure_generic(
        threshold_full,
        THRESHOLD_COLUMN_MAP,
        OUTPUTS["threshold_ever_full"]["svg"],
        OUTPUTS["threshold_ever_full"]["png"],
        note_text=THRESHOLD_TABLE_NOTE,
    )
    print()
    print(f"Saved: {OUTPUTS['threshold_ever_full']['svg']}")
    print(f"Saved: {OUTPUTS['threshold_ever_full']['png']}")

    threshold_total_bw = build_threshold_table(threshold_total_tables, POPULATION_MAP)
    save_table_figure_generic(
        threshold_total_bw,
        THRESHOLD_COLUMN_MAP,
        OUTPUTS["threshold_total_bw"]["svg"],
        OUTPUTS["threshold_total_bw"]["png"],
        note_text=THRESHOLD_TABLE_NOTE,
    )
    print()
    print(f"Saved: {OUTPUTS['threshold_total_bw']['svg']}")
    print(f"Saved: {OUTPUTS['threshold_total_bw']['png']}")

    threshold_total_full = build_threshold_table_full_pairs(threshold_total_tables)
    save_table_figure_generic(
        threshold_total_full,
        THRESHOLD_COLUMN_MAP,
        OUTPUTS["threshold_total_full"]["svg"],
        OUTPUTS["threshold_total_full"]["png"],
        note_text=THRESHOLD_TABLE_NOTE,
    )
    print()
    print(f"Saved: {OUTPUTS['threshold_total_full']['svg']}")
    print(f"Saved: {OUTPUTS['threshold_total_full']['png']}")

    probability_p_full = build_probability_table_full(probability_tables, FULL_POPULATION_MAP, "P")
    save_table_figure_generic(
        probability_p_full,
        MEASURE_COLUMN_MAP,
        OUTPUTS["probability_p_full"]["svg"],
        OUTPUTS["probability_p_full"]["png"],
        note_text=PROBABILITY_HIGHER_NOTE,
        label_width_override=1.35,
    )
    print()
    print(f"Saved: {OUTPUTS['probability_p_full']['svg']}")
    print(f"Saved: {OUTPUTS['probability_p_full']['png']}")

    probability_2pminus1_full = build_probability_table_full(probability_tables, FULL_POPULATION_MAP, "2Pminus1")
    save_table_figure_generic(
        probability_2pminus1_full,
        MEASURE_COLUMN_MAP,
        OUTPUTS["probability_2pminus1_full"]["svg"],
        OUTPUTS["probability_2pminus1_full"]["png"],
        note_text=PROBABILITY_DIRECTION_NOTE,
        label_width_override=1.35,
    )
    print()
    print(f"Saved: {OUTPUTS['probability_2pminus1_full']['svg']}")
    print(f"Saved: {OUTPUTS['probability_2pminus1_full']['png']}")

    figure1_plot_df, figure1_table_df = build_figure1_data(eligibility_heatmaps)
    save_figure1(
        figure1_plot_df,
        figure1_table_df,
        OUTPUTS["figure1"]["svg"],
        OUTPUTS["figure1"]["png"],
    )
    print()
    print(f"Saved: {OUTPUTS['figure1']['svg']}")
    print(f"Saved: {OUTPUTS['figure1']['png']}")


if __name__ == "__main__":
    main()
