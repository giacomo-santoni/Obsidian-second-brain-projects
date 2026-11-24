---
tags: [debug, dcache, pnfs, fermilab, scp]
date: 2026-11-24
---

# scp to /pnfs: "Permission denied"

## Error

```
scp: dest open "/pnfs/dune/persistent/users/gsantoni/beam_monitoring/data/fast_reco/nominal/CAF_fluxsyst.flat/0000000/sand_events_nominal_FHC.caffluxsyst.sand.full.sandstt1.0000020.CAF.flat.root": Permission denied
```

## Cause

`/pnfs/dune/...` is **dCache**, not a normal filesystem. dCache files are **write-once**: once a file is closed it can't be opened for writing again. `scp` tries to open the destination for writing and truncate it, and dCache refuses with a misleading `Permission denied`.

Most likely the destination file **already exists**. It can also be a 0-byte file left behind by an earlier `scp` that failed partway.

## Fixes

1. **If the file already exists, delete it first and copy again** (run on the gpvm):
   ```bash
   ls -l /pnfs/dune/persistent/users/gsantoni/.../0000000/
   rm /pnfs/dune/persistent/users/gsantoni/.../0000000/<file>.root
   ```

2. **Check that you own the directory.** If another user created it, you can't create files in it:
   ```bash
   ls -ld /pnfs/dune/persistent/users/gsantoni/.../0000000
   ```

3. **Use a dCache-aware tool instead of scp** (recommended):
   ```bash
   ifdh cp -D localfile.root /pnfs/dune/persistent/users/gsantoni/.../0000000/
   ```
   or `xrdcp` with a valid token (`htgettoken` first). `xrdcp -f` can overwrite a file.

   Workflow: `scp` the file to your home directory or `/exp/dune/data` on the gpvm, then `ifdh cp` it into `/pnfs`.

## Related
- [[Beam monitoring]]
