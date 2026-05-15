"""
Compatibility entry point for jthlab/bim.

The BIM 0.1.1 CLI converts scipy.optimize result arrays with ``float(out.x)``.
On current NumPy versions, one-element arrays must be converted explicitly.
This wrapper keeps the installed BIM package unmodified and patches only the
two BIM estimators used by the CLI before delegating to ``bim.cli``.
"""

def _scalarize_x(result):
    x = result.x
    if hasattr(x, "shape") and tuple(x.shape) == (1,):
        result.x = x.item()
    return result


def _patch_predict_methods():
    from bim.Bimbalance import bSFS, bTree

    if not getattr(bSFS.predict, "_sweep_scalar_patch", False):
        original_bsfs_predict = bSFS.predict

        def bsfs_predict(self, *args, **kwargs):
            return _scalarize_x(original_bsfs_predict(self, *args, **kwargs))

        bsfs_predict._sweep_scalar_patch = True
        bSFS.predict = bsfs_predict

    if not getattr(bTree.predict, "_sweep_scalar_patch", False):
        original_btree_predict = bTree.predict

        def btree_predict(self, *args, **kwargs):
            return _scalarize_x(original_btree_predict(self, *args, **kwargs))

        btree_predict._sweep_scalar_patch = True
        bTree.predict = btree_predict


def main():
    _patch_predict_methods()
    from bim.cli import cli

    cli(prog_name="bim")


if __name__ == "__main__":
    main()
