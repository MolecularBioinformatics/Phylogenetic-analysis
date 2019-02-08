TaxFinder
=========

Before the TaxFinder module can be used, a taxonomy database must be created. This is done by running `updateTaxonomy.sh`. The script will take up to an hour and will create about 6 GB data. During creation of the database, up to 20 GB disk space is needed.

In order to conveniently use TaxFinder as modules in scripts, create a symlink in some Python import path, e.g.:

```
ln -s /home/username/scripts/phylogenetics/taxfinder/taxfinder.py /home/username/.local/lib/pythonX.X/site-packages
```

Adjust the paths and the python version before running the above command. After this step, you can then import TaxFinder in your script and use its functions. Initializing TaxFinder in a script may take several seconds.

```python
from taxfinder import TaxFinder

TF = TaxFinder()

TF.getLineage(9606)
```
