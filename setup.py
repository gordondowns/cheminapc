#import setuptools
#
#with open("README.md", "r") as fh:
#    long_description = fh.read()
#
#setuptools.setup(
#    name="chemin-apc",
#    version="0.1",
#    author="Gordon Downs",
#    author_email="gdowns@email.arizona.edu",
#    description="X-ray Diffraction data analysis software, developed for the CheMin team of NASA's Mars Science Laboratory",
#    long_description=long_description,
#    long_description_content_type="text/markdown",
#    url="https://github.com/pypa/sampleproject",
#    packages=setuptools.find_packages(),
#    classifiers=[
#        "Programming Language :: Python :: 2",
#        "License :: OSI Approved :: GNU General Public License",
#        "Operating System :: OS Independent",
#    ],
#)

from distutils.core import setup
setup(
  name = 'cheminapc',         # How you named your package folder (MyLib)
  packages = ['cheminapc'],   # Chose the same as "name"
  version = '0.1',      # Start with a small number and increase it with every change you make
  license = 'gpl-3.0',        # Chose a license from here: https://help.github.com/articles/licensing-a-repository
  description = "X-ray Diffraction data analysis software, developed for the CheMin team of NASA's Mars Science Laboratory",   # Give a short description about your library
  author = 'Gordon Downs',                   # Type in your name
  author_email = 'gdowns@email.arizona.edu',      # Type in your E-Mail
  url = 'https://github.com/gordiedowns/cheminapc',   # Provide either the link to your github or to your website
  download_url = 'https://github.com/user/reponame/archive/v_01.tar.gz',    # I explain this later on
  keywords = ['XRD', 'CheMin'],   # Keywords that define your package best
  install_requires = [            # I get to this in a second
          'numpy',
          'scipy',
          'matplotlib',
      ],
  classifiers = [
    'Development Status :: 3 - Alpha',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Intended Audience :: Developers',      # Define that your audience are developers
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: GNU General Public License v3.0',   # Again, pick a license
    'Programming Language :: Python :: 2',      #Specify which pyhton versions that you want to support
    "Operating System :: OS Independent",
  ],
)