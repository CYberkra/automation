"""Prepare auditable pixel labels; never infer lithology, metres, or run gprMax.

The PNG and H5 outputs are private, provisional segmentation candidates, not
simulator inputs. Specify colours explicitly; unknown RGB values stay unknown
in exact labels and receive a flagged nearest-colour suggestion separately.
"""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
from PIL import Image


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(source, out, colours, border_colour=None):
    if out.exists():
        raise ValueError('Refuse to overwrite a previous preparation')
    source_hash = digest(source)
    with Image.open(source) as image:
        if image.mode != 'RGB':
            raise ValueError('Only opaque RGB images are supported')
        rgb = np.array(image)
    palette = np.asarray(colours, dtype=np.int32)
    if (palette.ndim != 2 or palette.shape[1] != 3 or len(palette) < 1 or
            len(palette) > 32767 or np.any(palette < 0) or np.any(palette > 255)
            or len(np.unique(palette, axis=0)) != len(palette)):
        raise ValueError('Palette must contain unique RGB triples in 0..255')
    original_shape = list(rgb.shape[:2])
    crop = None
    if border_colour is not None:
        border = np.all(rgb == np.asarray(border_colour), axis=-1)
        expected = np.zeros(rgb.shape[:2], dtype=bool)
        expected[-1, :] = True
        expected[:, -1] = True
        if not np.array_equal(border, expected):
            raise ValueError('Border colour is not confined to last row/column')
        # This crop is a declared candidate; source image remains unchanged.
        rgb = rgb[:-1, :-1]
        crop = dict(last_row=True, last_column=True, rgb=list(border_colour),
                    removed_pixels=int(expected.sum()), role='UNCONFIRMED_BORDER_CANDIDATE')
    distances = np.stack([
        np.sum((rgb.astype(np.int32)-colour)**2, axis=-1)
        for colour in palette
    ], axis=-1)
    best = distances.min(axis=-1)
    suggested = distances.argmin(axis=-1).astype(np.int16)
    exact = np.where(best == 0, suggested, -1).astype(np.int16)
    tied = (distances == best[..., None]).sum(axis=-1) > 1
    cleaned = palette[suggested].astype(np.uint8)
    if digest(source) != source_hash:
        raise ValueError('Source changed during preparation')
    out.mkdir(parents=True)
    labels = out/'image_labels_not_solver_geometry.h5'
    with h5py.File(labels, 'x') as h5:
        h5.create_dataset('exact_colour_labels', data=exact, compression='gzip')
        h5.create_dataset('nearest_colour_suggestion', data=suggested, compression='gzip')
        h5.create_dataset('changed_pixel_mask', data=best != 0, compression='gzip')
        h5.create_dataset('tied_nearest_colour_mask', data=tied, compression='gzip')
        h5.create_dataset('palette_rgb', data=palette.astype(np.uint8))
        h5.attrs['status'] = 'PIXEL_SEGMENTATION_DRAFT_NOT_GPRMAX_GEOMETRY'
        h5.attrs['axes'] = 'image_row_down,image_column_right'
        h5.attrs['units'] = 'pixels; physical scale not assigned'
        h5.attrs['source_sha256'] = source_hash
    with h5py.File(labels, 'r') as h5:
        assert np.array_equal(h5['exact_colour_labels'][:], exact)
        assert np.array_equal(h5['nearest_colour_suggestion'][:], suggested)
        assert np.array_equal(h5['changed_pixel_mask'][:], best != 0)
        assert np.array_equal(h5['tied_nearest_colour_mask'][:], tied)
    candidate = out/'palette_clean_candidate.png'
    Image.fromarray(cleaned).save(candidate)
    with Image.open(candidate) as saved:
        assert np.array_equal(np.array(saved), cleaned)
    result = dict(
        status='PREPARED_PIXEL_LABELS_AWAITING_LEGEND_SCALE_AND_BORDER_REVIEW',
        source_sha256=source_hash, script_sha256=digest(Path(__file__)),
        original_shape_row_column=original_shape,
        candidate_shape_row_column=list(exact.shape), border_candidate_crop=crop,
        unique_source_colours=int(len(np.unique(rgb.reshape(-1, 3), axis=0))),
        palette_rgb=palette.tolist(),
        exact_pixel_counts=[int((exact == i).sum()) for i in range(len(palette))],
        suggested_pixel_counts=[int((suggested == i).sum()) for i in range(len(palette))],
        unknown_exact_pixels=int((exact < 0).sum()),
        nearest_tie_pixels=int(tied.sum()),
        maximum_nearest_rgb_distance=float(np.sqrt(best.max())),
        material_names=None, physical_dimensions_m=None, air_colour=None,
        user_direction_note='User says lower part is ground/air; legend pending',
        orientation='Source image direction retained; no geological flip performed',
        official_geometry_axes='For future V4 export: XYZ, not legacy ZYX',
        outputs_sha256={p.name:digest(p) for p in (labels, candidate)},
        calls_solver=False, originals_unchanged=digest(source) == source_hash)
    (out/'summary.json').write_text(json.dumps(result, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2, ensure_ascii=False))


def triple(value):
    try:
        values = tuple(int(item) for item in value.split(','))
    except ValueError as exc:
        raise argparse.ArgumentTypeError('Expected R,G,B integers') from exc
    if len(values) != 3 or any(item < 0 or item > 255 for item in values):
        raise argparse.ArgumentTypeError('Expected R,G,B in 0..255')
    return values


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--image', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--colour', type=triple, action='append', required=True)
    parser.add_argument('--border-candidate', type=triple)
    args = parser.parse_args()
    prepare(args.image, args.out, args.colour, args.border_candidate)
