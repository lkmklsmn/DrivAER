import setuptools

setuptools.setup(
    name="DrivAER",
    version="0.0.1",
    description="Identify driving transcriptional regulators in single cell embeddings",
    author='Lukas Simon, Fangfang Yan',
    author_email="lkmklsmn@gmail.com",
    packages=['DrivAER'],
    install_requires=['sklearn','scanpy','anndata',
                      'pandas','seaborn','matplotlib','dca','tensorflow<=1.15.2'
                      ],
    package_data={'DrivAER': ['data/*.txt','annotations/*.gmt','annotations/*.tsv']}
)
