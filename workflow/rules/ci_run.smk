## CI run rules
## Mode B only: seed MACS2 peaks from the full run into each subsample's output directory
## so that the ABC peak-calling rules are skipped (ruleorder prefers these seeding rules).
##
## In mode A ("full"), no rules here are active — the main pipeline rules handle everything.

if MODE == "abc":

    rule seed_peaks_narrowpeak:
        ## Copy full run's narrowPeak file into the subsample output directory.
        ## Takes priority over abc_call_macs_peaks via ruleorder below.
        output:
            narrowpeak = os.path.join(SCRATCH_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak")
        params:
            src_dir = lambda wildcards: os.path.join(
                FULL_RUN_DIR, wildcards.biosample.rsplit("__sub", 1)[0], "Peaks"
            )
        run:
            import shutil
            os.makedirs(os.path.dirname(output.narrowpeak), exist_ok=True)
            # Copy all files in the Peaks directory (covers any auxiliary files ABC may need)
            for fname in os.listdir(params.src_dir):
                shutil.copy2(
                    os.path.join(params.src_dir, fname),
                    os.path.join(os.path.dirname(output.narrowpeak), fname)
                )

    rule seed_peaks_sorted:
        ## Provide the sorted narrowPeak from the full run.
        ## Takes priority over abc_sort_narrowpeaks via ruleorder below.
        input:
            os.path.join(SCRATCH_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak")
        output:
            sorted = os.path.join(SCRATCH_DIR, "{biosample}", "Peaks", "macs2_peaks.narrowPeak.sorted")
        params:
            src = lambda wildcards: os.path.join(
                FULL_RUN_DIR, wildcards.biosample.rsplit("__sub", 1)[0],
                "Peaks", "macs2_peaks.narrowPeak.sorted"
            )
        run:
            import shutil
            shutil.copy2(params.src, output.sorted)

    ruleorder: seed_peaks_narrowpeak > abc_call_macs_peaks
    ruleorder: seed_peaks_sorted > abc_sort_narrowpeaks
