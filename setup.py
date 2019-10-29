from setuptools import setup, find_packages

setup(
    name="genofunk",
    version="0.0.0",
    packages=find_packages(),
    url="https://github.com/rmcolq/genofunk",
    license="MIT",
    entry_points={"console_scripts": ["genofunk = genofunk.__main__:main"]},
    test_suite="nose.collector",
    tests_require=["nose >= 1.3"],
    install_requires=[
        "biopython>=1.70",
        "pandas>=0.24.2",
        "parasail>=1.1.17",
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3 :: Only",
        "License :: OSI Approved :: MIT License",
    ],
)
