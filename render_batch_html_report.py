#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import html
import os
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Sequence
from urllib.parse import quote


STATUS_ORDER = {
    "good_no_red_flags": 0,
    "mismatch_needs_review": 1,
    "analysis_failed": 2,
    "samplesheet_incomplete": 3,
    "samplesheet_extraction_failed": 4,
}

STATUS_LABELS = {
    "good_no_red_flags": "Good",
    "mismatch_needs_review": "Needs Review",
    "analysis_failed": "Analysis Failed",
    "samplesheet_incomplete": "Samplesheet Incomplete",
    "samplesheet_extraction_failed": "Extraction Failed",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render an HTML summary report from an IGVF batch TSV report."
    )
    parser.add_argument("--report-tsv", required=True, help="Batch TSV report to render")
    parser.add_argument("--output-html", help="Output HTML path; defaults to report.html next to the TSV")
    parser.add_argument("--title", default="IGVF Batch Validation Report", help="HTML report title")
    return parser.parse_args()


def load_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sort_rows(rows: Sequence[Dict[str, str]]) -> List[Dict[str, str]]:
    return sorted(
        rows,
        key=lambda row: (
            STATUS_ORDER.get(row.get("status", ""), 99),
            row.get("lab", ""),
            row.get("analysis_set_accession", ""),
        ),
    )


def compact_text(value: str, *, empty: str = "") -> str:
    token = str(value or "").strip()
    return token if token else empty


def rel_href(target: str, base_dir: Path) -> str:
    path = Path(target)
    if not path.exists():
        return ""
    relative = os.path.relpath(path.resolve(), base_dir.resolve())
    return quote(relative.replace(os.sep, "/"))


def render_link(label: str, target: str, base_dir: Path) -> str:
    href = rel_href(target, base_dir)
    if not href:
        return ""
    return f'<a href="{href}">{html.escape(label)}</a>'


def render_links(row: Dict[str, str], base_dir: Path) -> str:
    links = []
    for label, key in (
        ("index", "index_path"),
        ("summary", "analysis_summary_path"),
        ("samplesheet", "samplesheet_path"),
        ("stderr", "stderr_log"),
        ("stdout", "stdout_log"),
    ):
        target = compact_text(row.get(key, ""))
        if not target:
            continue
        link = render_link(label, target, base_dir)
        if link:
            links.append(link)
    return " | ".join(links) if links else ""


def render_details(summary: str, body: str, *, open_default: bool = False) -> str:
    safe_summary = html.escape(summary)
    safe_body = html.escape(body)
    if not body:
        return ""
    open_attr = " open" if open_default else ""
    return f"<details{open_attr}><summary>{safe_summary}</summary><div class=\"detail-body\">{safe_body}</div></details>"


def render_status_cards(counts: Counter[str]) -> str:
    cards = []
    for status in sorted(counts, key=lambda item: STATUS_ORDER.get(item, 99)):
        cards.append(
            "<div class=\"status-card\">"
            f"<div class=\"status-label {html.escape(status)}\">{html.escape(STATUS_LABELS.get(status, status))}</div>"
            f"<div class=\"status-count\">{counts[status]}</div>"
            "</div>"
        )
    return "".join(cards)


def render_lab_table(rows: Sequence[Dict[str, str]]) -> str:
    by_lab: Dict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        by_lab[compact_text(row.get("lab", ""), empty="Unknown")][row.get("status", "")] += 1
    body = []
    for lab in sorted(by_lab):
        status_counts = by_lab[lab]
        body.append(
            "<tr>"
            f"<td>{html.escape(lab)}</td>"
            f"<td>{sum(status_counts.values())}</td>"
            f"<td>{status_counts.get('good_no_red_flags', 0)}</td>"
            f"<td>{status_counts.get('mismatch_needs_review', 0)}</td>"
            f"<td>{status_counts.get('analysis_failed', 0)}</td>"
            f"<td>{status_counts.get('samplesheet_extraction_failed', 0)}</td>"
            "</tr>"
        )
    return "".join(body)


def render_rows(rows: Sequence[Dict[str, str]], base_dir: Path) -> str:
    body = []
    for row in rows:
        accession = compact_text(row.get("analysis_set_accession", ""))
        analysis_set_path = compact_text(row.get("analysis_set_path", ""))
        portal_link = ""
        if analysis_set_path:
            portal_link = f"https://api.data.igvf.org{analysis_set_path}"
        assay = compact_text(row.get("summary", ""), empty="NA")
        simplified = compact_text(row.get("simplified_sample_summary", ""))
        sample_summary = compact_text(row.get("sample_summary", ""))
        selection_mode = compact_text(row.get("selection_mode", ""), empty="NA")
        selection_note = compact_text(row.get("selection_note", ""))
        flags = compact_text(row.get("comparison_flags", ""), empty="NA")
        status = compact_text(row.get("status", ""), empty="NA")
        note = compact_text(row.get("status_note", ""), empty="NA")
        links = render_links(row, base_dir)

        accession_cell = html.escape(accession)
        if portal_link:
            accession_cell = f'<a href="{html.escape(portal_link)}">{html.escape(accession)}</a>'

        assay_html = html.escape(assay)
        if simplified:
            assay_html += render_details("sample", simplified)
        if sample_summary:
            assay_html += render_details("full sample summary", sample_summary)

        selection_html = html.escape(selection_mode)
        if selection_note:
            selection_html += render_details("selection note", selection_note)

        note_html = html.escape(note)
        if len(note) > 140:
            note_html = render_details(note[:140] + "...", note, open_default=False)

        body.append(
            "<tr>"
            f"<td>{accession_cell}</td>"
            f"<td>{html.escape(compact_text(row.get('lab', ''), empty='NA'))}</td>"
            f"<td><span class=\"status-pill {html.escape(status)}\">{html.escape(STATUS_LABELS.get(status, status))}</span></td>"
            f"<td>{assay_html}</td>"
            f"<td>{selection_html}</td>"
            f"<td>{html.escape(flags)}</td>"
            f"<td>{note_html}</td>"
            f"<td>{links}</td>"
            "</tr>"
        )
    return "".join(body)


def render_html(title: str, rows: Sequence[Dict[str, str]], source_tsv: Path, output_html: Path) -> str:
    counts = Counter(row.get("status", "") for row in rows)
    lab_table = render_lab_table(rows)
    table_rows = render_rows(rows, output_html.parent)
    generated_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{html.escape(title)}</title>
  <style>
    :root {{
      --bg: #efe9df;
      --panel: #fffaf4;
      --line: #decfbd;
      --text: #261d16;
      --muted: #6d6256;
      --accent: #8f4f19;
      --good: #dff2df;
      --good-text: #1e6a33;
      --review: #f5ead0;
      --review-text: #8b5f11;
      --fail: #f7dfda;
      --fail-text: #9a3426;
      --extract: #ece6f5;
      --extract-text: #5c4978;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      font-family: "Avenir Next", "Segoe UI", sans-serif;
      color: var(--text);
      background:
        radial-gradient(circle at top left, rgba(172, 121, 60, 0.16), transparent 28%),
        linear-gradient(180deg, #f4eee5 0%, var(--bg) 100%);
    }}
    main {{
      width: min(1480px, calc(100vw - 32px));
      margin: 24px auto 48px;
    }}
    .hero, .panel {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 24px;
      box-shadow: 0 18px 40px rgba(81, 53, 26, 0.08);
    }}
    .hero {{
      padding: 28px 30px;
    }}
    .hero p {{
      margin: 8px 0 0;
      color: var(--muted);
    }}
    .grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(170px, 1fr));
      gap: 14px;
      margin-top: 18px;
    }}
    .status-card {{
      background: #f8f2ea;
      border: 1px solid var(--line);
      border-radius: 18px;
      padding: 16px 18px;
    }}
    .status-label {{
      font-size: 12px;
      font-weight: 700;
      letter-spacing: 0.08em;
      text-transform: uppercase;
    }}
    .status-count {{
      margin-top: 10px;
      font-size: 32px;
      font-weight: 800;
    }}
    .panel {{
      margin-top: 18px;
      padding: 22px 24px;
    }}
    .toolbar {{
      display: flex;
      gap: 12px;
      flex-wrap: wrap;
      margin-top: 12px;
      align-items: center;
    }}
    .toolbar input {{
      width: min(360px, 100%);
      padding: 10px 12px;
      border-radius: 12px;
      border: 1px solid var(--line);
      background: #fff;
      font: inherit;
    }}
    .toolbar button {{
      padding: 10px 14px;
      border-radius: 999px;
      border: 1px solid var(--line);
      background: #fff7ed;
      color: var(--text);
      font: inherit;
      cursor: pointer;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      margin-top: 12px;
    }}
    th, td {{
      text-align: left;
      vertical-align: top;
      padding: 10px 12px;
      border-bottom: 1px solid #eadfce;
      font-size: 14px;
    }}
    th {{
      font-size: 12px;
      text-transform: uppercase;
      letter-spacing: 0.06em;
      color: var(--muted);
      position: sticky;
      top: 0;
      background: var(--panel);
    }}
    a {{
      color: #2b5ab3;
      text-decoration: none;
    }}
    .status-pill {{
      display: inline-flex;
      padding: 6px 10px;
      border-radius: 999px;
      font-size: 12px;
      font-weight: 700;
      text-transform: uppercase;
      white-space: nowrap;
    }}
    .good_no_red_flags {{
      background: var(--good);
      color: var(--good-text);
    }}
    .mismatch_needs_review {{
      background: var(--review);
      color: var(--review-text);
    }}
    .analysis_failed {{
      background: var(--fail);
      color: var(--fail-text);
    }}
    .samplesheet_incomplete,
    .samplesheet_extraction_failed {{
      background: var(--extract);
      color: var(--extract-text);
    }}
    details {{
      margin-top: 6px;
    }}
    summary {{
      cursor: pointer;
      color: var(--accent);
    }}
    .detail-body {{
      margin-top: 8px;
      white-space: pre-wrap;
      color: var(--muted);
    }}
    .table-wrap {{
      overflow: auto;
      max-height: 70vh;
      border: 1px solid var(--line);
      border-radius: 18px;
      background: #fffdf9;
    }}
    .meta {{
      margin-top: 10px;
      font-size: 13px;
      color: var(--muted);
    }}
  </style>
</head>
<body>
  <main>
    <section class="hero">
      <h1>{html.escape(title)}</h1>
      <p>Generated from <code>{html.escape(str(source_tsv))}</code> on {html.escape(generated_at)}.</p>
      <div class="grid">{render_status_cards(counts)}</div>
      <div class="meta">Rows: {len(rows)}. Good rows in this batch: {counts.get('good_no_red_flags', 0)}.</div>
    </section>

    <section class="panel">
      <h2>By Lab</h2>
      <div class="table-wrap">
        <table>
          <thead>
            <tr>
              <th>Lab</th>
              <th>Total</th>
              <th>Good</th>
              <th>Needs Review</th>
              <th>Analysis Failed</th>
              <th>Extraction Failed</th>
            </tr>
          </thead>
          <tbody>{lab_table}</tbody>
        </table>
      </div>
    </section>

    <section class="panel">
      <h2>Batch Rows</h2>
      <div class="toolbar">
        <input id="searchBox" type="search" placeholder="Filter by accession, lab, assay, note" />
        <button data-status="all">All</button>
        <button data-status="good_no_red_flags">Good</button>
        <button data-status="mismatch_needs_review">Needs Review</button>
        <button data-status="analysis_failed">Analysis Failed</button>
        <button data-status="samplesheet_extraction_failed">Extraction Failed</button>
      </div>
      <div class="table-wrap">
        <table id="reportTable">
          <thead>
            <tr>
              <th>Accession</th>
              <th>Lab</th>
              <th>Status</th>
              <th>Assay / Sample</th>
              <th>Selection</th>
              <th>Flags</th>
              <th>Note</th>
              <th>Artifacts</th>
            </tr>
          </thead>
          <tbody>{table_rows}</tbody>
        </table>
      </div>
    </section>
  </main>
  <script>
    const buttons = document.querySelectorAll('[data-status]');
    const rows = Array.from(document.querySelectorAll('#reportTable tbody tr'));
    const searchBox = document.getElementById('searchBox');
    let activeStatus = 'all';

    function applyFilters() {{
      const query = (searchBox.value || '').toLowerCase();
      rows.forEach((row) => {{
        const rowText = row.innerText.toLowerCase();
        const statusText = row.querySelector('.status-pill')?.classList[1] || '';
        const statusOk = activeStatus === 'all' || statusText === activeStatus;
        const queryOk = !query || rowText.includes(query);
        row.style.display = statusOk && queryOk ? '' : 'none';
      }});
    }}

    buttons.forEach((button) => {{
      button.addEventListener('click', () => {{
        activeStatus = button.dataset.status || 'all';
        applyFilters();
      }});
    }});
    searchBox.addEventListener('input', applyFilters);
  </script>
</body>
</html>
"""


def main() -> None:
    args = parse_args()
    report_tsv = Path(args.report_tsv).resolve()
    if args.output_html:
        output_html = Path(args.output_html).resolve()
    else:
        output_html = report_tsv.with_suffix(".html")
    rows = sort_rows(load_rows(report_tsv))
    output_html.parent.mkdir(parents=True, exist_ok=True)
    output_html.write_text(
        render_html(args.title, rows, report_tsv, output_html),
        encoding="utf-8",
    )
    print(f"Wrote HTML report to {output_html}")


if __name__ == "__main__":
    main()
