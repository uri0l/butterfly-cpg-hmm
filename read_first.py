# Re-import libraries after code execution state reset
from Bio import SeqIO

fasta_path = "GCA_933228805.2_ilHipSeme1.2_genomic.fna"

with open(fasta_path, "r") as fp:
	for record in SeqIO.parse(fp, "fasta"):
		output_path = f"chromosome_{record.description[58]}.fna"
		zw_records = []
		desc_lower = record.description.lower()
		if "chromosome: z" in desc_lower or "chromosome: w" in desc_lower or "chromosome: 28" in desc_lower:
			zw_records.append(record)
			SeqIO.write(zw_records, output_path, "fasta")
			for record in zw_records:
				print(f"Found: {record.id} | {record.description[58, 59]} | Length: {len(record.seq)}")


