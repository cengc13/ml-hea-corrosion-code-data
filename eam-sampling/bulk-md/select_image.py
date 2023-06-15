from amp.preprocess.image_selection import FurthestPointSampling
from amp.utilities import hash_images

images = 'bulk_md_nvt_1700K_al20.traj'
images = hash_images(images)
fps = FurthestPointSampling(images, k=20, encoder='gaussian')
chosen_images = fps.search(calculate_dev=True, save_traj='chosen_bulk_al20.traj')
