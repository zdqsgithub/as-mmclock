import GEOparse
import sys
try:
    print("Testing GEOparse download...")
    gse = GEOparse.get_GEO(geo="GSE80761", destdir="./")
    print(f"Success! {len(gse.gsms)} samples.")
except Exception as e:
    print(f"Error: {e}")
