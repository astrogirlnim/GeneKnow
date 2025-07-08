"""
Test script to debug report_sections issue.
"""
from graph import run_pipeline
import json

# Run the pipeline
result = run_pipeline(
    "../test_R1.fastq.gz",
    {
        "patient_data": {
            "age": 45,
            "sex": "F"
        }
    }
)

print("=" * 60)
print("PIPELINE RESULT KEYS:")
print("=" * 60)
print(list(result.keys()))

print("\n" + "=" * 60)
print("REPORT SECTIONS:")
print("=" * 60)
if "report_sections" in result:
    if result["report_sections"]:
        print(json.dumps(result["report_sections"], indent=2))
    else:
        print("report_sections is empty: {}")
else:
    print("report_sections NOT FOUND in result")

print("\n" + "=" * 60)
print("COMPLETED NODES:")
print("=" * 60)
print(result.get("completed_nodes", []))

print("\n" + "=" * 60)
print("PIPELINE STATUS:")
print("=" * 60)
print(result.get("pipeline_status", "unknown")) 