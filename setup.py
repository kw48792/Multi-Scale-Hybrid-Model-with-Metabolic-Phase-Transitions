"""
Setup script for the Multi-Scale Hybrid Model package.
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="multiscale-hybrid-model",
    version="1.0.0",
    author="[Author Names]",
    author_email="[author@email.com]",
    description="Multi-Scale Hybrid Modeling for Cell Culture with Metabolic Phase Transitions",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/kw48792/Multi-Scale-Hybrid-Model-with-Metabolic-Phase-Transitions",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=6.0",
            "black>=21.0",
            "flake8>=3.9",
        ],
    },
    include_package_data=True,
    zip_safe=False,
)
