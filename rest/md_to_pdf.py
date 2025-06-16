###Runs the pandoc shell command to convert .md to .pdfs and deletes the .md and intermediate files
import os


def report_to_pdf(path_to_md, path_to_pdf, sample_id):
    ''' Converts the markdown report to pdf using pandoc.'''
    os.system(f"pandoc {path_to_md} -o {path_to_pdf} --pdf-engine=tectonic --template=resources/template.tex")
    print(f"{sample_id} report saved at {path_to_pdf}")

def main(path_to_md, path_to_pdf, temp_dir = "tmpdir_hiv"):
    
    sample_id = path_to_md.split("/")[-1].replace(".md", "").replace("report_", "")
    report_to_pdf(path_to_md, path_to_pdf, sample_id)
    os.remove(path_to_md)  # Remove the .md file after conversion
    print(f"Removed {path_to_md}")
    for temp_file in os.listdir(temp_dir):
            if temp_file.startswith(f"{sample_id}"):
                os.remove(os.path.join(temp_dir, temp_file))
    print(f"Removed temporary files for {sample_id}")


if __name__ == "__main__":
    main(
        snakemake.input.report,
        snakmake.output.report,
    )
