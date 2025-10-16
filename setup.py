from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="gtgmm",
    version="2.0.0",
    author="Kevin Song, John Zhang, Lei Ye, Jianyi Zhang",
    author_email="kmsong@uab.edu",
    description="GeneTerrain Gaussian Mixture Models: Topological Data Analysis for Gene Networks",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/kevinmsong/gtGMM",
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
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "pandas>=1.3.0",
        "scikit-learn>=1.0.0",
        "scipy>=1.7.0",
        "matplotlib>=3.4.0",
        "networkx>=2.6.0",
        "requests>=2.26.0",
        "gseapy>=1.0.0",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=3.0.0",
            "black>=22.0.0",
            "flake8>=4.0.0",
            "mypy>=0.950",
            "jupyter>=1.0.0",
            "ipykernel>=6.0.0",
            "sphinx>=4.0.0",
            "sphinx_rtd_theme>=1.0.0",
        ],
    },
)

