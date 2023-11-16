import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="blind_rt60",
    version="0.1.0-a0",
    author="Asaf Zorea",
    author_email="zoreasaf@gmail.com",
    description="The BlindRT60 algorithm is used to estimate the reverberation time (RT60) "
                "of a room based on the recorded audio signals from microphones",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license='MIT',
    url="https://github.com/nuniz/blind_rt60",
    packages=setuptools.find_packages(exclude=["tests", "*.tests", "*.tests.*", "tests.*", "tests.*"]),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        'License :: OSI Approved :: MIT License',
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=[
        "scipy",
        "numpy",
        "matplotlib"
    ],
    extras_require={
        "dev": ["pyroomacoustics", "parameterized"],
    },
)
