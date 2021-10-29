import setuptools
 
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()
 
setuptools.setup(
    name="abageotests",
    version="0.0.5",
    author="WANG Hailin",
    author_email="hailin.wang@connect.polyu.hk",
    description="A package to generate and execute python scripts of abaqus models of geotechnical tests",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Hailin-Wang/abageotests",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)