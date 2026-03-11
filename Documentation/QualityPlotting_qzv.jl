### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ bccf3096-2775-4f7b-98ed-38f696df11fc
# cell-pkgs
begin
	import Pkg
	Pkg.add(["CairoMakie", "CSV", "DataFrames", "Statistics", "ZipFile", "JSON3", "PlutoUI"])
	using CairoMakie, CSV, DataFrames, Statistics, ZipFile, JSON3, PlutoUI
end

# ╔═╡ 2f1fe006-e64d-489a-a60e-3f5be5ce0b42
begin
	### A Pluto.jl notebook ###
	# v0.20.3
	
	using Markdown
	using InteractiveUtils
end

# ╔═╡ 36febd4d-0248-4073-9d34-eee78beb0abf
# cell-title
md"""
# QIIME2 .qzv QC Visualizer
Plot per-base quality, per-sequence quality, read counts, and denoising statistics
from your QIIME2 pipeline output files.
"""

# ╔═╡ 2233d87d-73f1-4cd8-a3a6-e610e4bf3479
# cell-paths-header
md"""
## Step 1 — Set Your File Paths
Edit `base_dir` and the file lists below to match your local paths.
"""

# ╔═╡ 4ce7b460-c0ae-4dea-9774-216be6a0b58a
# cell-paths-def
begin
	# base_dir = "/Users/jacksadler/Desktop/project_thesis/slurm_outputs/test_302"
	base_dir = "/Volumes/ROSALIND/cluster_outputs/full_12SV5"
	
	demux_files = [
		joinpath(base_dir, "demux.qzv"),
		joinpath(base_dir, "illumina_adapter_trimmed.qzv"),
		joinpath(base_dir, "adapter_and_primerp1.qzv"),
		joinpath(base_dir, "cutadapt_12SV5_finished.qzv"),
	]

	denoising_file = joinpath(base_dir, "denoising_stats.qzv")
	outdir         = joinpath(base_dir, "plots")
	mkpath(outdir)

	md"✅ Paths configured. Output directory: `$outdir`"
end

# ╔═╡ dcf73f32-0945-4da6-9b00-999d7155b306
# cell-utils
begin
	function extract_qzv(qzv_path::String)
		tmpdir = mktempdir()
		zf = ZipFile.Reader(qzv_path)
		for f in zf.files
			outpath = joinpath(tmpdir, f.name)
			mkpath(dirname(outpath))
			if !endswith(f.name, "/")
				write(outpath, read(f))
			end
		end
		close(zf)
		return tmpdir
	end

	function find_file(root::String, pattern::String)
		for (dirpath, _, files) in walkdir(root)
			for f in files
				occursin(pattern, f) && return joinpath(dirpath, f)
			end
		end
		return nothing
	end

	step_label(p::String) = join(uppercasefirst.(split(replace(basename(p), r"\.qzv$" => ""), "_")), " ")
	step_stem(p::String)  = replace(basename(p), r"\.qzv$" => "")

	md"✅ Utility functions loaded"
end

# ╔═╡ 6b271135-7726-41bd-b81b-086ce9b0db23
# cell-parsers
begin
	function parse_perbase(tmpdir::String)
		fwd_path = find_file(tmpdir, "forward-seven-number-summaries.tsv")
		rev_path = find_file(tmpdir, "reverse-seven-number-summaries.tsv")
		fwd = isnothing(fwd_path) ? nothing : CSV.read(fwd_path, DataFrame; delim='\t')
		rev = isnothing(rev_path) ? nothing : CSV.read(rev_path, DataFrame; delim='\t')
		return fwd, rev
	end

	function parse_counts(tmpdir::String)
		path = find_file(tmpdir, "per-sample-fastq-counts.tsv")
		isnothing(path) && return nothing
		return CSV.read(path, DataFrame; delim='\t')
	end

	function parse_denoising(tmpdir::String)
		path = find_file(tmpdir, "metadata.tsv")
		isnothing(path) && return nothing
		lines      = readlines(path)
		data_lines = filter(l -> !startswith(l, "#"), lines)
		isempty(data_lines) && return nothing
		return CSV.read(IOBuffer(join(data_lines, "\n")), DataFrame; delim='\t', header=1)
	end

	md"✅ Parsing functions loaded"
end


# ╔═╡ 3c710e28-0e18-4faf-8907-3e7fe9f009bc
# cell-plotfns
begin
	function plot_perbase_quality(fwd, rev, label, stem, outdir)
		fig = Figure(size = (1000, 500))
		ax  = Axis(fig[1,1],
				   title  = "$label — Per-Base Quality Scores",
				   xlabel = "Position in Read (bp)",
				   ylabel = "Phred Quality Score",
				   yticks = 0:5:45)

		hspan!(ax, 28, 45; color = (:green,  0.07))
		hspan!(ax, 20, 28; color = (:orange, 0.07))
		hspan!(ax,  0, 20; color = (:red,    0.07))

		function extract_pb(df)
			positions  = 1:length(names(df)[2:end])
			row_labels = lowercase.(string.(df[:, 1]))
			get_row(tag) = begin
				idx = findfirst(r -> occursin(tag, r), row_labels)
				isnothing(idx) ? fill(NaN, length(positions)) : Float64.(Vector(df[idx, 2:end]))
			end
			return positions, get_row("50"), get_row("25"), get_row("75")
		end

		pos_f, med_f, q25_f, q75_f = extract_pb(fwd)
		band!(ax, pos_f, q25_f, q75_f; color = (:steelblue, 0.25))
		lines!(ax, pos_f, med_f; color = :steelblue, linewidth = 2, label = "Forward (R1)")

		if !isnothing(rev)
			pos_r, med_r, q25_r, q75_r = extract_pb(rev)
			band!(ax, pos_r, q25_r, q75_r; color = (:tomato, 0.25))
			lines!(ax, pos_r, med_r; color = :tomato, linewidth = 2, label = "Reverse (R2)")
		end

		axislegend(ax; position = :rb)
		ylims!(ax, 0, 45)
		save(joinpath(outdir, "$(stem)_perbase_quality.pdf"), fig)
		return fig
	end

	function plot_perseq_quality(fwd, rev, label, stem, outdir)
		fig = Figure(size = (900, 450))
		ax  = Axis(fig[1,1],
				   title  = "$label — Per-Sequence Median Quality",
				   xlabel = "Median Phred Quality Score",
				   ylabel = "Count of Positions")

		function get_median_row(df)
			idx = findfirst(r -> occursin("50", r), lowercase.(string.(df[:, 1])))
			isnothing(idx) && return nothing
			return Float64.(Vector(df[idx, 2:end]))
		end

		m1 = get_median_row(fwd)
		!isnothing(m1) && hist!(ax, m1; bins = 0:1:42, color = (:steelblue, 0.6), label = "Forward (R1)")

		if !isnothing(rev)
			m2 = get_median_row(rev)
			!isnothing(m2) && hist!(ax, m2; bins = 0:1:42, color = (:tomato, 0.6), label = "Reverse (R2)")
		end

		axislegend(ax; position = :lt)
		save(joinpath(outdir, "$(stem)_perseq_quality.pdf"), fig)
		return fig
	end

	function plot_read_counts(counts, label, stem, outdir)
		fig = Figure(size = (max(600, nrow(counts) * 18), 500))
		ax  = Axis(fig[1,1],
				   title               = "$label — Read Counts per Sample",
				   xlabel              = "Sample",
				   ylabel              = "Read Count",
				   xticklabelrotation = π/3,
				   xticklabelsize     = 9)

		col_names = names(counts)
		fwd_col   = findfirst(c -> occursin("forward", lowercase(c)) || occursin("r1", lowercase(c)), col_names)
		rev_col   = findfirst(c -> occursin("reverse", lowercase(c)) || occursin("r2", lowercase(c)), col_names)
		samples   = string.(counts[:, first(col_names)])
		xs        = 1:nrow(counts)

		!isnothing(fwd_col) && barplot!(ax, xs, Float64.(counts[:, fwd_col]);
			color = (:steelblue, 0.75), label = "Forward (R1)", dodge = fill(1, nrow(counts)))
		!isnothing(rev_col) && barplot!(ax, xs, Float64.(counts[:, rev_col]);
			color = (:tomato, 0.75), label = "Reverse (R2)", dodge = fill(2, nrow(counts)))

		ax.xticks = (xs, samples)
		axislegend(ax; position = :rt)
		save(joinpath(outdir, "$(stem)_read_counts.pdf"), fig)
		return fig
	end

	function plot_denoising_stats(df, stem, outdir)
		fig = Figure(size = (1000, 550))
		ax  = Axis(fig[1,1],
				   title               = "DADA2 Denoising Statistics",
				   xlabel              = "Sample",
				   ylabel              = "Read Count",
				   xticklabelrotation = π/3,
				   xticklabelsize     = 9)

		col_names = lowercase.(names(df))
		samples   = string.(df[:, names(df)[1]])
		xs        = 1:nrow(df)

		stage_map = [
			("input",        "Input",        :gray60),
			("filtered",     "Filtered",     :steelblue),
			("denoised",     "Denoised",     :mediumseagreen),
			("merged",       "Merged",       :mediumpurple),
			("non-chimeric", "Non-Chimeric", :tomato),
		]

		for (i, (key, lbl, col)) in enumerate(stage_map)
			m = findfirst(c -> occursin(key, c), col_names)
			isnothing(m) && continue
			barplot!(ax, xs, Float64.(df[:, names(df)[m]]);
					 dodge = fill(i, nrow(df)), color = (col, 0.75), label = lbl)
		end

		ax.xticks = (xs, samples)
		axislegend(ax; position = :rt)
		save(joinpath(outdir, "$(stem)_denoising_stats.pdf"), fig)
		return fig
	end

	md"✅ Plotting functions loaded"
end

# ╔═╡ 76446cfc-ce9c-4aa6-bedf-8d6e594d7eb0
# cell-dropdown-header
md"""
## Step 2 — Select Trimming Step to Visualize
"""


# ╔═╡ 31188a3a-3c33-42ff-94c1-427fa797f593
# cell-dropdown
@bind selected_demux Select([f => step_label(f) for f in demux_files])


# ╔═╡ c5587ca4-4907-4cdd-82eb-f8aa46462b7f
# cell-selected-label
md"Selected: **$(step_label(selected_demux))**"

# ╔═╡ 10a83130-2a2e-49aa-9fc1-147a9079ee1e
let
    tmpdir = extract_qzv(demux_files[1])
    for (dirpath, dirs, files) in walkdir(tmpdir)
        for f in files
            println(joinpath(dirpath, f))
        end
    end
    rm(tmpdir; recursive = true)
end

# ╔═╡ 3966ac35-5ed3-49bd-9ac9-185d724ab4f6
# cell-perbase-header
md"## Step 3 — Per-Base Quality"

# ╔═╡ 93190aaa-66b0-4842-bffc-d73380d08dd3
# cell-perbase
let
	label  = step_label(selected_demux)
	stem   = step_stem(selected_demux)
	tmpdir = extract_qzv(selected_demux)
	fwd, rev = parse_perbase(tmpdir)
	fig = !isnothing(fwd) ? plot_perbase_quality(fwd, rev, label, stem, outdir) :
		  md"⚠️ Per-base quality data not found."
	rm(tmpdir; recursive = true)
	fig
end

# ╔═╡ b999fd9e-8f49-4b19-951e-d32416646e3b
# cell-perseq-header
md"## Step 4 — Per-Sequence Quality"


# ╔═╡ ca5f577c-6892-41c4-a9ec-cf6836d2988d
# cell-perseq
let
	label  = step_label(selected_demux)
	stem   = step_stem(selected_demux)
	tmpdir = extract_qzv(selected_demux)
	fwd, rev = parse_perbase(tmpdir)
	fig = !isnothing(fwd) ? plot_perseq_quality(fwd, rev, label, stem, outdir) :
		  md"⚠️ Per-sequence quality data not found."
	rm(tmpdir; recursive = true)
	fig
end

# ╔═╡ 11b378d4-f7d5-499f-a305-68752892f249
# cell-counts-header
md"## Step 5 — Read Counts per Sample"

# ╔═╡ 137479b1-1312-402e-b2b0-f2efb087d68f
# cell-counts
let
	label  = step_label(selected_demux)
	stem   = step_stem(selected_demux)
	tmpdir = extract_qzv(selected_demux)
	counts = parse_counts(tmpdir)
	fig = !isnothing(counts) ? plot_read_counts(counts, label, stem, outdir) :
		  md"⚠️ Read count data not found."
	rm(tmpdir; recursive = true)
	fig
end

# ╔═╡ b1c34d31-417f-49e6-be03-6e3021540d2b
# cell-gc-header
md"## Step 6 — GC Content"

# ╔═╡ ddbe427e-5ff1-4764-b0d4-c18445a3e060
# cell-gc
md"""
!!! warning "GC Content Not Available in .qzv"
    QIIME2's `demux summarize` does not export per-read GC content to its data directory.
    Run **FastQC** or **MultiQC** on your raw `.fastq.gz` files for GC content plots,
    using the `plot_fastqc.jl` script.
"""


# ╔═╡ ece282f5-4b95-4dd9-82be-33ba688af237
# cell-denoising-header
md"## Step 7 — Denoising Statistics"

# ╔═╡ ef5ef9b1-ff81-4c4b-8c13-f8aae175c843
let
    tmpdir = extract_qzv(denoising_file)
    for (dirpath, _, files) in walkdir(tmpdir)
        for f in files
            println(joinpath(dirpath, f))
        end
    end
    rm(tmpdir; recursive = true)
end

# ╔═╡ f1e7ba50-002c-493a-9603-9f6ee2746c3b
# cell-denoising
let
	stem   = step_stem(denoising_file)
	tmpdir = extract_qzv(denoising_file)
	df     = parse_denoising(tmpdir)
	fig = !isnothing(df) ? plot_denoising_stats(df, stem, outdir) :
		  md"⚠️ Could not parse denoising stats."
	rm(tmpdir; recursive = true)
	fig
end

# ╔═╡ c84dc5ce-1c9b-476a-95aa-f25ba6248a4d
# cell-export-header
md"""
## Step 8 — Export All Steps to PDF
Run this cell to generate PDFs for **all** demux files at once.
"""


# ╔═╡ 97b30448-98a2-455f-bb6c-756d50a043ee
# cell-export
begin
	for qzv in demux_files
		label  = step_label(qzv)
		stem   = step_stem(qzv)
		tmpdir = extract_qzv(qzv)
		fwd, rev = parse_perbase(tmpdir)
		counts   = parse_counts(tmpdir)
		!isnothing(fwd)    && plot_perbase_quality(fwd, rev, label, stem, outdir)
		!isnothing(fwd)    && plot_perseq_quality(fwd, rev, label, stem, outdir)
		!isnothing(counts) && plot_read_counts(counts, label, stem, outdir)
		rm(tmpdir; recursive = true)
	end

	tmpdir = extract_qzv(denoising_file)
	df     = parse_denoising(tmpdir)
	!isnothing(df) && plot_denoising_stats(df, step_stem(denoising_file), outdir)
	rm(tmpdir; recursive = true)

	md"✅ All PDFs exported to `$outdir`"
end

# ╔═╡ Cell order:
# ╟─2f1fe006-e64d-489a-a60e-3f5be5ce0b42
# ╟─bccf3096-2775-4f7b-98ed-38f696df11fc
# ╟─36febd4d-0248-4073-9d34-eee78beb0abf
# ╟─2233d87d-73f1-4cd8-a3a6-e610e4bf3479
# ╠═4ce7b460-c0ae-4dea-9774-216be6a0b58a
# ╠═dcf73f32-0945-4da6-9b00-999d7155b306
# ╠═6b271135-7726-41bd-b81b-086ce9b0db23
# ╠═3c710e28-0e18-4faf-8907-3e7fe9f009bc
# ╟─76446cfc-ce9c-4aa6-bedf-8d6e594d7eb0
# ╠═31188a3a-3c33-42ff-94c1-427fa797f593
# ╠═c5587ca4-4907-4cdd-82eb-f8aa46462b7f
# ╠═10a83130-2a2e-49aa-9fc1-147a9079ee1e
# ╟─3966ac35-5ed3-49bd-9ac9-185d724ab4f6
# ╟─93190aaa-66b0-4842-bffc-d73380d08dd3
# ╟─b999fd9e-8f49-4b19-951e-d32416646e3b
# ╟─ca5f577c-6892-41c4-a9ec-cf6836d2988d
# ╟─11b378d4-f7d5-499f-a305-68752892f249
# ╟─137479b1-1312-402e-b2b0-f2efb087d68f
# ╟─b1c34d31-417f-49e6-be03-6e3021540d2b
# ╟─ddbe427e-5ff1-4764-b0d4-c18445a3e060
# ╟─ece282f5-4b95-4dd9-82be-33ba688af237
# ╟─ef5ef9b1-ff81-4c4b-8c13-f8aae175c843
# ╟─f1e7ba50-002c-493a-9603-9f6ee2746c3b
# ╟─c84dc5ce-1c9b-476a-95aa-f25ba6248a4d
# ╟─97b30448-98a2-455f-bb6c-756d50a043ee
