# Code Upload Guide

This guide helps you add your journal paper code to this repository.

## Steps to Upload Your Code

### 1. Add Source Code

Place your Python code files in the `src/` directory:

```bash
src/
├── __init__.py                 # Create this if it doesn't exist
├── your_main_model.py          # Your main model implementation
├── preprocessing.py            # Data preprocessing code
├── training.py                 # Model training code
└── ... (other modules)
```

### 2. Add Data Files

Place your data files in the `data/` directory:
- Raw experimental data → `data/raw/`
- Processed data → `data/processed/`

**Note**: For large files (>100 MB), consider:
- Using Git LFS
- Hosting on Zenodo/figshare and adding download links

### 3. Add Example Scripts

Create example notebooks or scripts in `examples/`:
- `examples/01_data_preprocessing.ipynb`
- `examples/02_model_training.ipynb`
- `examples/03_results_visualization.ipynb`

### 4. Update Documentation

#### 4.1. Update Main README

Edit `README.md` to include:
- Specific installation instructions
- Usage examples
- Data description
- Citation information (once paper is published)

#### 4.2. Update Requirements

Edit `requirements.txt` with exact versions of packages you used:
```txt
numpy==1.21.5
pandas==1.3.4
# etc.
```

#### 4.3. Update Citation File

Edit `CITATION.cff` with:
- Author names and ORCID IDs
- Paper title and DOI (once published)
- Publication year

#### 4.4. Update License

Edit `LICENSE` file if you want to use a different license:
- MIT (currently set) - Permissive, good for academic code
- Apache 2.0 - Similar to MIT with patent grant
- GPL-3.0 - Copyleft license
- CC-BY-4.0 - For data and documentation

### 5. Add Additional Documentation

In `docs/`, create:
- `model_description.md` - Mathematical formulation
- `reproduction_guide.md` - How to reproduce paper results
- Any additional technical documentation

### 6. Test Your Code

Before committing:
1. Test installation: `pip install -r requirements.txt`
2. Run examples to ensure they work
3. Check that all file paths are relative to the repository root

### 7. Organize by Paper Sections (Optional)

If your paper has multiple experiments, consider organizing:

```
examples/
├── figure_1_simulation/
├── figure_2_validation/
├── table_1_comparison/
└── supplementary_analysis/
```

## Best Practices

### Code Quality
- Add docstrings to functions and classes
- Include comments for complex algorithms
- Use meaningful variable names
- Follow PEP 8 style guide

### Reproducibility
- Document random seeds used
- Include version numbers for all dependencies
- Provide example data or data generation scripts
- Document computational environment (OS, hardware)

### Data Management
- Include data dictionaries (variable descriptions)
- Document units for all measurements
- Provide sample data if full dataset is proprietary
- Use standard file formats (CSV, HDF5, etc.)

### Version Control
- Make small, logical commits
- Write clear commit messages
- Don't commit large binary files directly
- Use `.gitignore` to exclude temporary files

## Checklist Before Submission

- [ ] All source code is in `src/`
- [ ] Data files are organized in `data/`
- [ ] At least one working example is provided
- [ ] `requirements.txt` is complete and tested
- [ ] README is updated with specific project information
- [ ] CITATION.cff has author information
- [ ] LICENSE is appropriate for your needs
- [ ] Code runs without errors
- [ ] Examples reproduce key paper results
- [ ] Documentation explains how to use the code

## Questions?

If you need help:
1. Check existing issues in the repository
2. Create a new issue with your question
3. Contact the repository maintainer

## Example Commit Workflow

```bash
# Add your files
git add src/your_model.py
git add data/processed/experiment_data.csv
git add examples/demo.ipynb

# Commit with a clear message
git commit -m "Add hybrid model implementation and training data"

# Push to GitHub
git push origin main
```

## After Upload

Once your code is uploaded:
1. Test the repository by cloning it fresh and following your own README
2. Ask a colleague to try reproducing results
3. Update the paper manuscript with the GitHub URL
4. Consider creating a release when the paper is accepted
5. Add a DOI using Zenodo integration

---

**Remember**: The goal is to make your research reproducible and accessible to other researchers in your field!
