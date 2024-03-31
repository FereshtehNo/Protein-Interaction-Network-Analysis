# Load required packages
library(bio3d)
library(igraph)

# Step1:  Set the project directory
project_dir <- "D:/Phd-classes/Final-project"
input_dir <- file.path(project_dir, "ProcessedData")
output_dir <- file.path(project_dir, "AnalysisResults")

# Step2: Create input and output directories if they don't exist
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Step3: Define the function to save specific lines and save files
process_and_save <- function(input_file, output_file, start_line, end_line) {
  lines_desired_range <- scan(input_file, what = character(), skip = start_line - 1, nlines = end_line - start_line + 1, sep = "\n")
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  
  # Save in AnalysisResults directory
  writeLines(lines_desired_range, con = output_file)
  
  # Save in ProcessedData directory
  processed_output_file <- file.path(input_dir, basename(output_file))
  dir.create(dirname(processed_output_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines_desired_range, con = processed_output_file)
}



# Step4: Define the function to calculate distance matrix for all atoms
calculate_distance_matrix <- function(pdb_file_name, output_file_name) {
  pdb_input_path <- file.path(input_dir, pdb_file_name)
  pdb_output_path <- file.path(output_dir, pdb_file_name)
  
  pdb_lines <- readLines(pdb_input_path)
  
  # Extract the atomic coordinates from: https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
  coords <- matrix(nrow = length(pdb_lines), ncol = 3)
  row_count <- 0
  
  for (line in pdb_lines) {
    record_type <- substr(line, 1, 6)
    if (record_type == "ATOM  " || record_type == "HETATM") {
      x <- as.numeric(substr(line, 31, 38))
      y <- as.numeric(substr(line, 39, 46))
      z <- as.numeric(substr(line, 47, 54))
      row_count <- row_count + 1
      coords[row_count, ] <- c(x, y, z)
    }
  }
  
  # Remove empty rows
  coords <- coords[complete.cases(coords), ]
  
  # Calculate the distance matrix
  dist_matrix <- as.matrix(dist(coords))
  
  # Create a new column with row numbers
  row_numbers <- matrix(seq_len(nrow(dist_matrix)))
  
  # Combine the row numbers with the distance matrix
  dist_matrix_with_row_numbers <- cbind(row_numbers, dist_matrix)
  
  # Write the distance matrix to a text file in the input and output directories
  input_output_file_name <- paste0(output_file_name, ".txt")
  output_path_input <- file.path(input_dir, input_output_file_name)
  output_path_output <- file.path(output_dir, input_output_file_name)
  
  write.table(dist_matrix_with_row_numbers, file = output_path_input, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(dist_matrix_with_row_numbers, file = output_path_output, sep = "\t", quote = FALSE, row.names = FALSE)
}

# Step5: Define the function to calculate distance matrix for carbon alpha

calculate_ca_distance_matrix <- function(pdb_file_name, output_file_name) {
  pdb_input_path <- file.path(input_dir, pdb_file_name)
  pdb_output_path <- file.path(output_dir, pdb_file_name)
  
  pdb_lines <- readLines(pdb_input_path)
  
  # Extract the atomic coordinates for Carbon alpha (CA) atoms
  ca_coords <- matrix(nrow = length(pdb_lines), ncol = 3)
  row_count <- 0
  
  for (line in pdb_lines) {
    record_type <- substr(line, 1, 6)
    atom_name <- trimws(substr(line, 13, 16))  # Remove leading/trailing spaces
    if (record_type == "ATOM  " && atom_name == "CA") {
      x <- as.numeric(substr(line, 31, 38))
      y <- as.numeric(substr(line, 39, 46))
      z <- as.numeric(substr(line, 47, 54))
      row_count <- row_count + 1
      ca_coords[row_count, ] <- c(x, y, z)
    }
  }
  
  ca_coords <- ca_coords[complete.cases(ca_coords), ]
  
  # Calculate the distance matrix for Carbon alpha (CA) atoms
  dist_matrix <- as.matrix(dist(ca_coords))
  
  # Write the distance matrix to a text file in the input and output directories
  input_output_file_name <- paste0(output_file_name, ".txt")
  output_path_input <- file.path(input_dir, input_output_file_name)
  output_path_output <- file.path(output_dir, input_output_file_name)
  
  write.table(dist_matrix, file = output_path_input, sep = "\t", quote = FALSE)
  write.table(dist_matrix, file = output_path_output, sep = "\t", quote = FALSE)
}


# Step6: Define the function to create network for all atoms and carbon alpha
create_network <- function(input_file, output_file, cutoff_distance = 6) {
  input_path <- file.path(input_dir, input_file)
  
  dist_matrix <- as.matrix(read.table(input_path))
  
  n_atoms <- nrow(dist_matrix)
  network <- dist_matrix
  
  for (i in 1:n_atoms) {
    for (j in 1:n_atoms) {
      if (dist_matrix[i, j] > cutoff_distance) {
        network[i, j] <- 0
        network[j, i] <- 0   # The matrix is symmetric for an undirected network
      }
    }
  }
  
  # Save the network matrix to both ProcessedData and AnalysisResults directories
  output_path_input <- file.path(input_dir, output_file)
  output_path_output <- file.path(output_dir, output_file)
  
  write.table(network, file = output_path_input, sep = "\t", quote = FALSE)
  write.table(network, file = output_path_output, sep = "\t", quote = FALSE)
}


# Step7: Define the function to create and plot network
create_and_plot_network <- function(input_file) {
  input_path <- file.path(input_dir, input_file)
  
  network <- as.matrix(read.table(input_path))
  
  # Create an igraph graph from the adjacency matrix
  g <- graph.adjacency(network, mode = "undirected", weighted = TRUE)
  
  # Round the edge weights to a specific number of decimal places
  rounded_weights <- round(E(g)$weight, digits = 2)
  
  # Generate the output PDF filename based on the input filename
  pdf_file <- gsub("network_matrix_", "", input_file)
  pdf_file <- gsub(".txt", "_network_plot.pdf", pdf_file)
  
  # Plot the graph using Fruchterman-Reingold layout and save as PDF
  pdf(file.path(output_dir, pdf_file))
  plot(g, layout = layout_with_fr, edge.label = rounded_weights)
  dev.off()  # Close the PDF device
}

# Download PDB files
pdb_ids <- c("4YDF", "1RPI", "1LCS")

for (pdb_id in pdb_ids) {
  pdb_url <- paste0("https://files.rcsb.org/download/", pdb_id, ".pdb")
  download.file(pdb_url, file.path(input_dir, paste0(pdb_id, ".pdb")), mode = "wb")
}


#Save specific domains
# 1RPI cases
process_and_save(
  file.path(input_dir, "1RPI.pdb"),
  file.path(output_dir, "HIV-DOMAIN-A-ALL-ATOMS.pdb"),
  start_line = 333,
  end_line = 1115
)

process_and_save(
  file.path(input_dir, "1RPI.pdb"),
  file.path(output_dir, "HIV-DOMAIN-A-ACTIVE-SITE.pdb"),
  start_line = 536,
  end_line = 554
)

process_and_save(
  file.path(input_dir, "1RPI.pdb"),
  file.path(output_dir, "HIV-DOMAIN-A-FLAP.pdb"),
  start_line = 703,
  end_line = 778
)

process_and_save(
  file.path(input_dir, "1RPI.pdb"),
  file.path(output_dir, "HIV-DOMAIN-A-LOOP.pdb"),
  start_line = 960,
  end_line = 1000
)

# 4YDF cases
process_and_save(
  file.path(input_dir, "4YDF.pdb"),
  file.path(output_dir, "HTLV-1-DOMAIN-A-ALL-ATOMS.pdb"),
  start_line = 407,
  end_line = 1292
)

process_and_save(
  file.path(input_dir, "4YDF.pdb"),
  file.path(output_dir, "HTLV-1-DOMAIN-A-ACTIVE-SITE.pdb"),
  start_line = 638,
  end_line = 699
)

process_and_save(
  file.path(input_dir, "4YDF.pdb"),
  file.path(output_dir, "HTLV-1-DOMAIN-A-FLAP.pdb"),
  start_line = 819,
  end_line = 886
)

process_and_save(
  file.path(input_dir, "4YDF.pdb"),
  file.path(output_dir, "HTLV-1-DOMAIN-A-LOOP.pdb"),
  start_line = 1117,
  end_line = 1180
)

#create distance matrix for all atoms
calculate_distance_matrix("HIV-DOMAIN-A-ALL-ATOMS.pdb", "distance_matrix_HIV-DOMAIN-A-ALL-ATOMS")
calculate_distance_matrix("HTLV-1-DOMAIN-A-ALL-ATOMS.pdb", "distance_matrix_HTLV-1-DOMAIN-A-ALL-ATOMS")
calculate_distance_matrix("HIV-DOMAIN-A-ACTIVE-SITE.pdb", "distance_matrix_HIV-DOMAIN-A-ACTIVE-SITE")
calculate_distance_matrix("HTLV-1-DOMAIN-A-ACTIVE-SITE.pdb", "distance_matrix_HTLV-1-DOMAIN-A-ACTIVE-SITE")
calculate_distance_matrix("HIV-DOMAIN-A-FLAP.pdb", "distance_matrix_HIV-DOMAIN-A-FLAP")
calculate_distance_matrix("HTLV-1-DOMAIN-A-FLAP.pdb", "distance_matrix_HTLV-1-DOMAIN-A-FLAP")
calculate_distance_matrix("HIV-DOMAIN-A-LOOP.pdb", "distance_matrix_HIV-DOMAIN-A-LOOP")
calculate_distance_matrix("HTLV-1-DOMAIN-A-LOOP.pdb", "distance_matrix_HTLV-1-DOMAIN-A-LOOP")

#create distance matrix for carbon alpha
calculate_ca_distance_matrix("HIV-DOMAIN-A-ALL-ATOMS.pdb", "ca_distance_matrix_HIV-DOMAIN-A-ALL-ATOMS")
calculate_ca_distance_matrix("HTLV-1-DOMAIN-A-ALL-ATOMS.pdb", "ca_distance_matrix_HTLV-1-DOMAIN-A-ALL-ATOMS")
calculate_ca_distance_matrix("HIV-DOMAIN-A-ACTIVE-SITE.pdb", "ca_distance_matrix_HIV-DOMAIN-A-ACTIVE-SITE")
calculate_ca_distance_matrix("HTLV-1-DOMAIN-A-ACTIVE-SITE.pdb", "ca_distance_matrix_HTLV-1-DOMAIN-A-ACTIVE-SITE")
calculate_ca_distance_matrix("HIV-DOMAIN-A-FLAP.pdb", "ca_distance_matrix_HIV-DOMAIN-A-FLAP")
calculate_ca_distance_matrix("HTLV-1-DOMAIN-A-FLAP.pdb", "ca_distance_matrix_HTLV-1-DOMAIN-A-FLAP")
calculate_ca_distance_matrix("HIV-DOMAIN-A-LOOP.pdb", "ca_distance_matrix_HIV-DOMAIN-A-LOOP")
calculate_ca_distance_matrix("HTLV-1-DOMAIN-A-LOOP.pdb", "ca_distance_matrix_HTLV-1-DOMAIN-A-LOOP")

# create network for all atoms
create_network("distance_matrix_HIV-DOMAIN-A-ALL-ATOMS.txt", "network_matrix_HIV-DOMAIN-A-ALL-ATOMS.txt")
create_network("distance_matrix_HTLV-1-DOMAIN-A-ALL-ATOMS.txt", "network_matrix_HTLV-1-DOMAIN-A-ALL-ATOMS.txt")
create_network("distance_matrix_HIV-DOMAIN-A-ACTIVE-SITE.txt", "network_matrix_HIV-DOMAIN-A-ACTIVE-SITE.txt")
create_network("distance_matrix_HTLV-1-DOMAIN-A-ACTIVE-SITE.txt", "network_matrix_HTLV-1-DOMAIN-A-ACTIVE-SITE.txt")
create_network("distance_matrix_HIV-DOMAIN-A-FLAP.txt", "network_matrix_HIV-DOMAIN-A-FLAP.txt")
create_network("distance_matrix_HTLV-1-DOMAIN-A-FLAP.txt", "network_matrix_HTLV-1-DOMAIN-A-FLAP.txt")
create_network("distance_matrix_HIV-DOMAIN-A-Loop.txt", "network_matrix_HIV-DOMAIN-A-Loop.txt")
create_network("distance_matrix_HTLV-1-DOMAIN-A-Loop.txt", "network_matrix_HTLV-1-DOMAIN-A-Loop.txt")


# create network for carbon alpha
create_network("ca_distance_matrix_HIV-DOMAIN-A-ALL-ATOMS.txt", "network_matrix_ca_distance_matrix_HIV-DOMAIN-A-ALL-ATOMS.txt")
create_network("ca_distance_matrix_HTLV-1-DOMAIN-A-ALL-ATOMS.txt", "network_matrix_ca_distance_matrix_HTLV-1-DOMAIN-A-ALL-ATOMS.txt")
create_network("ca_distance_matrix_HIV-DOMAIN-A-ACTIVE-SITE.txt", "network_matrix_ca_distance_matrix_HIV-DOMAIN-A-ACTIVE-SITE.txt")
create_network("ca_distance_matrix_HTLV-1-DOMAIN-A-ACTIVE-SITE.txt", "network_matrix_ca_distance_matrix_HTLV-1-DOMAIN-A-ACTIVE-SITE.txt")
create_network("ca_distance_matrix_HIV-DOMAIN-A-FLAP.txt", "network_matrix_ca_distance_matrix_HIV-DOMAIN-A-FLAP.txt")
create_network("ca_distance_matrix_HTLV-1-DOMAIN-A-FLAP.txt", "network_matrix_ca_distance_matrix_HTLV-1-DOMAIN-A-FLAP.txt")
create_network("ca_distance_matrix_HIV-DOMAIN-A-Loop.txt", "network_matrix_ca_distance_matrix_HIV-DOMAIN-A-Loop.txt")
create_network("ca_distance_matrix_HTLV-1-DOMAIN-A-Loop.txt", "network_matrix_ca_distance_matrix_HTLV-1-DOMAIN-A-Loop.txt")

# List of input file names to draw plots
input_files <- c(
  "network_matrix_ca_distance_matrix_HIV-DOMAIN-A-ALL-ATOMS.txt",
  "network_matrix_ca_distance_matrix_HTLV-1-DOMAIN-A-ALL-ATOMS.txt",
  "network_matrix_ca_distance_matrix_HIV-DOMAIN-A-ACTIVE-SITE.txt",
  "network_matrix_ca_distance_matrix_HTLV-1-DOMAIN-A-ACTIVE-SITE.txt",
  "network_matrix_ca_distance_matrix_HIV-DOMAIN-A-FLAP.txt",
  "network_matrix_ca_distance_matrix_HTLV-1-DOMAIN-A-FLAP.txt",
  "network_matrix_ca_distance_matrix_HIV-DOMAIN-A-Loop.txt",
  "network_matrix_ca_distance_matrix_HTLV-1-DOMAIN-A-Loop.txt",
  "network_matrix_HIV-DOMAIN-A-ALL-ATOMS.txt",
  "network_matrix_HTLV-1-DOMAIN-A-ALL-ATOMS.txt",
  "network_matrix_HIV-DOMAIN-A-ACTIVE-SITE.txt",
  "network_matrix_HTLV-1-DOMAIN-A-ACTIVE-SITE.txt",
  "network_matrix_HIV-DOMAIN-A-FLAP.txt",
  "network_matrix_HTLV-1-DOMAIN-A-FLAP.txt",
  "network_matrix_HIV-DOMAIN-A-Loop.txt",
  "network_matrix_HTLV-1-DOMAIN-A-Loop.txt"
)


for (input_file in input_files) {
  create_and_plot_network(input_file)
}

