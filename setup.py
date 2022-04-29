from setuptools import setup
import io

# read the contents of your README file
from os import path

this_directory = path.abspath(path.dirname(__file__))
with io.open(path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="rss",
    packages=["rss"],
    version="1.0.0",
    description="Radical Stability Score",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Shree Sowndarya SV, Peter St. John, Robert Paton",
    author_email="patonlab@colostate.edu",
    url="https://github.com/patonlab/RSS",
    keywords=["compchem", "Radical Stability", "informatics"],
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
    ],
    install_requires=["numpy", "pandas", "cclib", "dbstep"],
    python_requires=">=3.6",
    include_package_data=True,
    package_data={"": ["*.csv"]},
)
