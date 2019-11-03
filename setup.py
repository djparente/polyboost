import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="polyboost", # Replace with your own username
    version="1.0.2",
    author="Daniel J. Parente",
    author_email="dparente@kumc.edu",
    description="An enhanced genomic variant classifier",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/djparente/polyboost",
    packages=setuptools.find_packages(),
    install_requires=[
        'numpy',
        'xgboost'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
    include_package_data=True
)