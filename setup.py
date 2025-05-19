from setuptools import setup, find_packages

setup(
    name="dna_sequencing",
    version="0.1.0",
    description="DNA Sequencing Algorithms Toolkit",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Your Name",
    author_email="your.email@example.com",
    url="https://github.com/andrewbrowne3/Algorithms_for_DNA_sequencing",
    packages=find_packages(),
    install_requires=[
        "matplotlib>=3.4.0",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.6",
) 