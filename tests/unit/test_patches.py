from snakemake.caching.local import OutputFileCache as LocalOutputFileCache

from showyourwork import patches


def test_patch_snakemake_cache():
    # Check that the patch changes the method function changed
    old_store = LocalOutputFileCache.store
    old_fetch = LocalOutputFileCache.fetch
    patches.patch_snakemake_cache(None, None)
    assert old_store == LocalOutputFileCache.store
    assert old_fetch == LocalOutputFileCache.fetch
