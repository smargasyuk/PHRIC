import filecmp


rule test:
    input:
        pipeline_out = ["results/test_hg19/S16.tsv"],
        ref_out = ["resources/test_results/test_hg19/S16.tsv"]
    run:
        for f1, f2 in zip(input["pipeline_out"], input["ref_out"]):
            if not filecmp.cmp(f1, f2, shallow = False):
                print(f"Files {f1} and {f2} are different")
                raise ValueError
        print("Test passed")


rule clean_test:
    shell: """
rm -rf results/test_hg19/
"""