import setuptools

verstr = "0.1.0"

setuptools.setup(
    name="modwaveforms",
    version=verstr,
    author="Jose María Ezquiaga",
    author_email="jose.ezquiaga@nbi.ku.dk",
    description="Collection of modified waveforms for parameter estimation",
    packages=[
        "modwaveforms",
        "modwaveforms.utils",
    ],
    install_requires=[
        "numpy",
        "scipy",
        "mpmath",
        "bilby",
        "bilby_pipe",
        "lalsuite",
    ],
    python_requires='>=3.7',
)
