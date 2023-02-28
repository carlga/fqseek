from setuptools import setup
from fqseek.__init__ import __version__, __author__, __email__

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name="fqseek",
    version=__version__,
    license="MIT",
    author=__author__,
    author_email=__email__,
    maintainer=__author__,
    maintainer_email=__email__,
    description=(
        "fqseek: an integrated framework for reproducible and efficient processing "
        "of high-throughput sequencing data."
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=["fqseek"],
    package_data={"" : ["fqseek/*"]},
    include_package_data=True,
    zip_safe=False,
    python_requires=">=3.6",
    install_requires=[],
    # setup_requires=open("requirements.txt").read().strip().split("\n"),
    url="https://github.com/carlga/fqseek",
    keywords=["fqseek", "NGS", "high-throughput sequencing", "QC", "fastq"],
    entry_points={
        "console_scripts": ["fqseek=fqseek.main:main"],
    },
    classifiers=[
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)