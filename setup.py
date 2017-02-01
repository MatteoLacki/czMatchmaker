from distutils.core import setup

setup(
    name          = 'czMatchmaker',
    packages      = ['czMatchmaker'], # this must be the same as the name above
    version       = '0.1',
    description   = 'Maximum likelihood estimates of electron transfer reaction parameters.',
    author        = 'Mateusz Krzysztof Lacki',
    author_email  = 'matteo.lacki@gmail.com',
    url           = 'https://github.com/MatteoLacki/linearCounter',
    download_url  = 'https://github.com/MatteoLacki/linearCounter.git',
    keywords      = ['Matching c and z ions', 'Mass Spectrometry'],
    license       = 'GNU AFFERO GENERAL PUBLIC LICENSE v3',
    classifiers   = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry'
    ],
    install_requires = ['dplython' , 'pandas']
)
