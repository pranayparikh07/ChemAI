# 📦 Data Files Setup Guide

The large CSV data files are **compressed** to fit GitHub's file size limits.

## ⚡ First Time Setup (After Cloning)

Run this to decompress the data:

```bash
python decompress_data.py
```

This will automatically extract all `.csv.gz` files to their original uncompressed versions.

---

## 📊 Large Files

These files are stored as compressed `.csv.gz`:

| File | Size (Uncompressed) | Size (Compressed) | Location |
|------|-------------------|------------------|----------|
| bioactivity_edges.csv | 253.0 MB | 31.5 MB | TEAM_SHREYA/data/ |
| molecules_kg.csv | 4.3 MB | 0.7 MB | TEAM_SHREYA/data/ |
| proteins_kg.csv | 0.9 MB | 0.2 MB | TEAM_SHREYA/data/ |
| similarity_edges.csv | 0.1 MB | 0.0 MB | TEAM_SHREYA/data/ |
| clean_candidates.csv | 7.3 MB | 2.0 MB | data/ |

---

## 🗜️ Compression Ratio

Files are compressed with gzip level 9 (maximum):
- **bioactivity_edges.csv**: 88% reduction ✅
- **molecules_kg.csv**: 84% reduction ✅
- **proteins_kg.csv**: 81% reduction ✅
- **similarity_edges.csv**: 85% reduction ✅
- **clean_candidates.csv**: 73% reduction ✅

**Total saved**: ~232 MB (allowing safe GitHub storage)

---

## 🔄 How It Works

1. **In Repository**: Files stored as `.csv.gz` (compressed)
2. **After Cloning**: Run `python decompress_data.py`
3. **For Development**: Use uncompressed `.csv` files normally
4. **Before Pushing**: Files stay as `.csv.gz` automatically

---

## 📝 For Team Members

### After Cloning
```bash
git clone https://github.com/pranayparikh07/ChemAI.git
cd ChemAI
python decompress_data.py  # ← Run this!
```

### Running ChemAI
```bash
# Now use as normal - decompressed files are ready
python run_chemai.py
python TEAM_SHREYA/load_to_neo4j.py
```

### Before Committing
- Compressed files are already in `.gitignore` as `*.csv`
- Only `.csv.gz` files are tracked by Git
- Changes to uncompressed files won't be committed

---

## 🔧 Recompressing Data (If Needed)

If you modify the data files and want to recompress:

```python
python -c "
import gzip
from pathlib import Path

file_path = Path('TEAM_SHREYA/data/bioactivity_edges.csv')
output = file_path.with_suffix('.csv.gz')

with open(file_path, 'rb') as f_in:
    with gzip.open(output, 'wb', compresslevel=9) as f_out:
        f_out.write(f_in.read())
        
print(f'✓ Recompressed: {output.name}')
"
```

---

## 💾 Alternative: Regenerate Data

Instead of using compressed files, you can **regenerate the data** from ChEMBL:

```bash
cd TEAM_SHREYA
python load_to_neo4j.py
```

See `TEAM_SHREYA/` for data generation scripts.

---

## ❓ Troubleshooting

**Q: Decompression fails?**
```bash
# Verify .gz files exist
dir *.csv.gz /s
python decompress_data.py -v
```

**Q: Modified data is too large?**
- Recompress before pushing
- Or regenerate from source

**Q: Lost uncompressed files?**
```bash
# Re-decompress from .gz
python decompress_data.py
```

---

## 📖 See Also

- [DATA_SETUP.md](./DATA_SETUP.md) - This file
- [decompress_data.py](./decompress_data.py) - Decompression script
- [TEAM_SHREYA/](../TEAM_SHREYA/) - Data generation scripts
- [PROJECT_DOCUMENTATION_HUB/](../PROJECT_DOCUMENTATION_HUB/) - Full documentation

---

**Status**: ✅ All data files successfully compressed and ready for GitHub!
