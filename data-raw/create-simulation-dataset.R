# Create a toy dataset to demonstrate usage of the Magellan package.
setwd("~/git/Magellan")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(Seurat)
library(Matrix)
library(RANN)
library(splatterBatch) ## https://github.com/jordansquair/splatter_batch
library(scater)
library(Magellan)

# set the seed
set.seed(1000)

# set the number of spatial transcriptomes
## here, simulate just two sections for speed
n_slices = 2

# simulate spatial coordinates for 100 barcodes
n_rows = 30
n_cols = 10
coords = matrix(1, nrow = n_rows, ncol = n_cols) %>%
  as.data.frame() %>%
  mutate(X = row_number()) %>%
  gather(Y, value, -X) %>%
  mutate(Y = gsub("V", "", Y)) %>%
  type_convert() %>%
  # now, re-scale to put this in the same coordinate space/sizing as the
  # main spatial object
  mutate(Y = scales::rescale(Y, c(25, 336)),
         X = scales::rescale(X, c(27, 482))) %>%
  dplyr::rename(imagecol = X, imagerow = Y) %>%
  mutate(tissue = 1) %>%
  mutate(row = imagerow, col = imagecol) %>%
  dplyr::select(tissue, row, col, imagerow, imagecol) %>%
  mutate(group = as.numeric(cut(imagecol, n_slices))) %>%
  # jitter the points
  mutate(imagecol = jitter(imagecol, amount = 1),
         imagerow = jitter(imagerow, amount = 2)) %>%
  # arrange by group to make sure we get the cell types right
  arrange(group)

# now, simulate transcriptomes for each barcode
batch_sizes = table(coords$group) %>% as.numeric()
# load simulation parameters
params = readRDS("data-raw/parameters.rds")
# define perturbation response
de_faclocs = c(1, 5)
cell_type_labels = paste0("CellType", LETTERS[seq_len(n_slices)])
cell_list = list()
for (index in seq_len(n_slices)) {
  cells = splatterBatch::splatSimulateGroups(
    params = params,
    batchCells = batch_sizes[index],
    # keep the seed the same, to generate the same 'cell type' over and over.
    seed = 1,
    de.prob = 0.1,
    de.facLoc = de_faclocs[index],
    group.prob = c(0.5, 0.5), verbose = F
  ) %>% logNormCounts() %>% as.Seurat
  # randomly shuffle the order of these cells
  set.seed(index)
  meta = cells@meta.data %>% extract(sample(rownames(.)),)
  expr = GetAssayData(cells, slot = 'counts') %>% extract(, rownames(meta))
  cells = CreateSeuratObject(expr, meta.data = meta) %>% NormalizeData()
  cells$cell_type = cell_type_labels[index]
  cells$de_facLoc = de_faclocs[index]
  cell_list[[index]] = cells
}

# reset seed
set.seed(1000)

# merge sections into a single object
sc = merge(cell_list[[1]], cell_list[-1])
sc$orig.ident = 1
final_sc = CreateSeuratObject(counts = GetAssayData(sc), assay = 'Spatial')

# create coordinates
reference = readRDS(file.path(Sys.getenv("PROJECTS"), "spatial-master.rds"))
reference = readRDS(file.path(Sys.getenv("PROJECTS"), "sc-walk/spatial/seurat/sc-walk-spatial.rds"))
reference = readRDS("~/git/sc-walk/data/spatial/seurat/sc-walk-spatial.rds")
final_sc@images$slice1 = reference@images$slice1
final_sc@images$slice1@coordinates = coords
final_sc@images$slice1@image = as.matrix(1)
rownames(coords) = colnames(final_sc)

# add label to the meta data
final_sc$label = sc$Group
final_sc$intensity = sc$Group
coords$label = sc$Group
coords %<>% dplyr::rename(intensity = group)

# save
# saveRDS(final_sc, "data/demo.rds")
demo = final_sc
save(demo, file = "data/demo.rda")

# plot for readme
p1 = coords %>%
  ggplot(aes(x = imagerow, y = imagecol, fill = intensity)) +
  geom_point(shape = 21, size = 0.8, stroke = 0) +
  coord_fixed() +
  scale_fill_distiller(palette = 'RdGy',
                       name = 'Perturbation\nresponse',
                       breaks = c(1.1, 1.9),
                       limits = c(1, 2),
                       labels = c(1, 5)) +
  guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
  umap_theme() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'right',
        legend.justification = 'bottom',
        legend.key.width = unit(0.2, 'lines'),
        legend.key.height = unit(0.18, 'lines'))
p1

p2 = coords %>%
  ggplot(aes(x = imagerow, y = imagecol, fill = label)) +
  geom_point(shape = 21, size = 0.8, stroke = 0) +
  coord_fixed() +
  scale_fill_brewer(palette = 'Paired',
                       name = 'Label') +
  # guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
  umap_theme() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'right',
        legend.justification = 'bottom',
        legend.key.width = unit(0.2, 'lines'),
        legend.key.height = unit(0.18, 'lines'))

ground_truth = p1 | p2

# run magellan
results = navigate_space(demo, n_threads = 16)
saveRDS(results, "results.rds")
results = readRDS("results.rds")
auc = results$AUC$AUC
range = range(auc$auc)
lims = c(range[1] + 0.1 * diff(range),
         range[2] - 0.1 * diff(range))
p3 = cbind(coords, auc) %>%
  ggplot(aes(x = imagerow, y = imagecol, fill = auc)) +
  geom_point(shape = 21, size = 0.8, stroke = 0) +
  scale_fill_distiller(palette = 'RdGy',
                       name = 'AUC',
                       breaks = lims,
                       limits = range,
                       labels = round(range, digits = 2)) +
  coord_fixed() +
  guides(fill = guide_colorbar(frame.colour = 'black', ticks = FALSE)) +
  umap_theme() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'right',
        legend.justification = 'bottom',
        legend.key.width = unit(0.2, 'lines'),
        legend.key.height = unit(0.18, 'lines'),
        )
p3

ggsave("fig/toy-ground-truth.png", ground_truth, width = 10, height = 3, units = 'cm', dpi = 600)
ggsave("fig/toy-magellan.png", p3, width = 10, height = 3, units = 'cm', dpi = 600)
