# Dashboard development and review

This directory contains a static Quarto dashboard for the lymphoma projections.
It runs entirely in the reader's browser and does not require Shiny, a database,
or access to the fitted model object.

## Refresh the public dashboard data

From the repository root:

```bash
Rscript dashboard/export_dashboard_data.R
```

The export script reads the canonical aggregate/modelled CSVs in `output/`, runs
release checks, and writes:

- `dashboard/data/lymphoma-dashboard-data.csv`
- `dashboard/data/data-dictionary.csv`

These are generated public-facing artefacts. Do not edit them by hand.

## Preview locally

Install Quarto, then run this from the repository root:

```bash
quarto preview dashboard
```

Quarto prints a local address, normally `http://localhost:xxxx/`. Open it in a
browser. Changes to `index.qmd` or `styles.css` are rebuilt while preview is
running. Stop the preview with `Ctrl+C`.

To create a non-live build instead:

```bash
quarto render dashboard
```

The rendered site is written to `dashboard/_site/` and is intentionally ignored
by Git.

## Review checklist

1. Confirm each selector presents only valid combinations.
2. Compare the 2021 and 2045 values against the current manuscript tables.
3. Check the incidence/prevalence boundary and uncertainty labelling.
4. Review the NHL hierarchy wording for possible double-counting ambiguity.
5. Test the CSV and data-dictionary downloads.
6. Review on a narrow/mobile window and at 200% browser zoom.
7. Confirm that no restricted, patient-level, confidential, or unpublished
   source fields appear in the exported CSV.
8. Replace the draft status, citation and DOI only after the corresponding
   research outputs are final.

## Publish with GitHub Pages

The workflow at `.github/workflows/publish-dashboard.yml` builds and deploys the
site after a push to `main`. Before the first deployment, a repository
administrator must select **GitHub Actions** under **Settings → Pages → Build and
deployment → Source**.
