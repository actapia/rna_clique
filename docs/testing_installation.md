# Testing the installation

After you have downloaded RNA-clique and installed its dependencies, you can
test your installation using a script included in this repository.

First, make sure your rna-clique conda environment is active.

```bash
conda activate rna-clique
```

Then, if you are in the root the repository, you can run the following command
to begin the test script.

```bash
bash tests/verify_install/test_install.sh && echo "Success!"
```

The script generates a small test dataset and runs RNA-clique on the generated
data. On a modern desktop with one thread, the test should take around one
minute to complete. On machines with multiple threads, the test script should
take advantage of parallelism to complete the test more quickly.

If you ran the script with the above command and see "Success!," then the
installation was succesful. Otherwise, you will need to investigate the output
of the test script to see what failed and why.

If the test script fails despite having a correct installation, you should
submit a bug report on GitHub at
[https://github.com/actapia/rna_clique/issues](https://github.com/actapia/rna_clique/issues)
.

