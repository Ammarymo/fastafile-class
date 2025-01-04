from setuptools import setup, find_packages

setup(
    name="fastafile",
    version="1.0.0",
    author="Ammar Y. Mohamed",
    author_email="amar.add655@gmail.com",
    description="A Python class for handling FASTA files.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/Ammarymo/fastafile-class",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
