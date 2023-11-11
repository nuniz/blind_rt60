import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="blind_rt60",
    version="0.1.0-alpha",
    author="Asaf Zorea",
    author_email="zoreasaf@gmail.com",
    description="Python implementation of blind estimation of reverberation time (rt60)",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license='MIT',
    url="https://github.com/nuniz/pyBlindRT",
    packages=setuptools.find_packages(exclude=["tests", "*.tests", "*.tests.*", "tests.*", "tests.*"]),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        'License :: OSI Approved :: MIT License',
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=[
        "matplotlib",
        "numpy",
        "soundfile",
    ],
)