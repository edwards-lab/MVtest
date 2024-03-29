from setuptools import setup, find_packages, find_namespace_packages

version = {}
with open("meanvar/version.py") as fp:
    exec(fp.read(), version)

import os


# Use the README as the long description
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name="mvtest",
    version=version["__version__"],
    author="Eric Torstenson",
    author_email="eric.s.torstenson@vanderbilt.edu",
    url="https://github.com/edwards-lab/MVtest",
    download_url="https://github.com/edwards-lab/MVtest/archive/v1.5.0.tar.gz",
    packages=["meanvar", "meanvar.tests"],
    license="GPL",
    scripts=["mvtest.py", "scripts/mvmany.py"],
    description=["GWAS analysis"],
    install_requires=["scipy", "numpy", "libgwas"],
    long_description=read("README"),
    keywords=["GWAS", "genetic analysis"],
    test_suite="meanvar.tests",
    package_data={"meanvar/tests/bedfiles/": ["*"], "doc": ["*"]},
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Utilities",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Software Development :: Libraries",
        "Programming Language :: Python :: 3.7",
    ],
)
