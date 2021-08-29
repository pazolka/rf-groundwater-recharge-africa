import numpy as np
import rasterio.crs
import rasterio.mask
from rasterio.enums import Resampling
import fiona
from rasterio.io import MemoryFile
from affine import Affine
from pyproj import Proj, transform

def clipAndResample(raster, shape, out_file, scaling_factor, crs="epsg:4326", allTouched=False):
    '''
    Clip raster (.asc or .tif) to a .shp polygon
    Rescale and/or resample the clipped raster
    Save to an output file

    parameters:
    - raster: path to the raster
    - shape: path to the .shp polygon
    - out_file: path to the output raster
    - scaling_factor: scaling factor
    - crs: the input projection (required for an .asc file)
    - allTouched: keep pixels that intersect with the polygon
    '''
    # clip 
    out_image, out_meta = clip(raster, shape, crs=crs, allTouched=allTouched)

    # create in-memory cropped raster
    with MemoryFile() as memfile:
        with memfile.open(**out_meta) as dataset:  # Open as DatasetWriter
            dataset.write(out_image)
        with memfile.open() as dataset:  # Reopen as DatasetReader
            # resample data to target resolution
            data = dataset.read(
                out_shape=(
                    dataset.count,
                    int(dataset.height * scaling_factor),
                    int(dataset.width * scaling_factor)
                ),
                resampling=Resampling.bilinear
            )

            # scale image transform
            transform = dataset.transform * dataset.transform.scale(
                (dataset.width / data.shape[-1]),
                (dataset.height / data.shape[-2])
            )
            # final raster config
            profile = dataset.profile
            profile.update(transform=transform, driver='GTiff', 
                height=int(dataset.height * scaling_factor), width=int(dataset.width * scaling_factor))

    # write file
    with rasterio.open(out_file, "w", **profile) as dest:
        dest.write(data)


def clip(raster, shape, crs="epsg:4326", allTouched=False):
    '''
    Clip raster (.asc or .tif) to a .shp polygon
    
    parameters:
    - raster: path to the raster
    - shape: path to the .shp polygon
    - crs: the input projection (required for an .asc file)
    - allTouched: keep pixels that intersect with the polygon
    '''
    # read continent outline
    with fiona.open(shape, "r") as shapefile:
        shapes = [feature["geometry"] for feature in shapefile]

    with rasterio.open(raster, mode='r+') as src:
        # in case it's an .asc file, add crs manually
        if src.crs == None:
            src.crs = rasterio.crs.CRS({"init": crs})
        # mask out the ocean
        out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True, all_touched=allTouched)
        out_meta = src.meta

        # config
        out_meta.update({"height": out_image.shape[1], "width": out_image.shape[2],"transform": out_transform})
    return out_image, out_meta


def clipAndSave(raster, shape, out_file, crs="epsg:4326", allTouched=False):
    '''
    Clip raster (.asc or .tif) to a .shp polygon
    Save to an output file

    parameters:
    - raster: path to the raster
    - shape: path to the .shp polygon
    - out_file: path to the output raster
    - crs: the input projection (required for an .asc file)
    - allTouched: keep pixels that intersect with the polygon
    '''
    # clip
    out_image, out_meta = clip(raster, shape, crs=crs, allTouched=allTouched)
    # save
    with rasterio.open(out_file, "w", **out_meta) as dest:
        dest.write(out_image)
        
def getAllCoordinates(raster, grid=0.5):
    '''
    Get all coordinates in a given raster
    https://gis.stackexchange.com/a/129857

    returns:
       - coords: a list of all coordinates (lat, long)
       - T0: Affine transformation of that raster
       - shape: initial raster shape
    '''
    # Read raster
    with rasterio.open(raster, 'r') as src:
        T0 = src.transform  # upper-left pixel corner affine transform
        p1 = Proj(src.crs)
        A = src.read()  # pixel values

    # All rows and columns
    cols, rows = np.meshgrid(np.arange(A.shape[2]), np.arange(A.shape[1]))

    # Get affine transform for pixel centres
    T1 = T0 * Affine.translation(grid, grid)
    # Function to convert pixel row/column index (from 0) to easting/northing at centre
    rc2en = lambda r, c: (c, r) * T1
    
    # All eastings and northings (there is probably a faster way to do this)
    eastings, northings = np.vectorize(rc2en, otypes=[float, float])(rows, cols)
    
    # Project all longitudes, latitudes
    p2 = Proj(proj='latlong',datum='WGS84')
    longs, lats = transform(p1, p2, eastings, northings)
    # to list
    coords = [item for item in zip(lats.flat, longs.flat)]
    return coords, T0, lats.shape