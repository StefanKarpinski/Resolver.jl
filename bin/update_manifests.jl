#!/usr/bin/env julia
#
# Refresh bin/'s per-Julia manifests against the current registries.
#
#     julia bin/update_manifests.jl            refresh every manifest
#     julia bin/update_manifests.jl --check    report gaps, change nothing
#
# bin/ pins its dependencies in one manifest per supported Julia minor, and one
# of those pins -- HistoricalStdlibVersions -- is a snapshot of what each Julia
# release bundles. Pins do not rot loudly: they keep resolving, they just start
# describing a Julia that no longer exists (issue #103). bin/check_stdlibs.jl
# catches that after the fact; this script is how you fix it, and how you keep
# it from happening.
#
# Each manifest has to be written by the Julia it is for, since it records that
# Julia's own stdlib versions -- so this drives `julia +X.Y` through juliaup,
# once per supported minor. Channels that are not installed are reported rather
# than installed: adding a Julia is your call, not a script's.
#
# The supported minors are every released Julia matching bin/Project.toml's
# `[compat] julia`, taken from the version list Registries.jl already tracks. A
# minor that has no manifest yet gets one -- that is the other half of staying
# current, and the half nothing else would notice, since a Julia with no
# manifest of its own does not fail, it silently falls back to another's.
#
# `--check` reports that gap without filling it, which is what CI runs: the
# check needs a current list of Julia releases, not a current set of pins, so it
# needs neither juliaup nor the other Julias to be installed.

push!(empty!(LOAD_PATH), @__DIR__)

import Pkg

include("Registries.jl")

const BIN = @__DIR__

# A minor is carried as a VersionNumber so it can be matched against compat,
# but it names a series, not a release: print it as one.
series(minor::VersionNumber) = "$(minor.major).$(minor.minor)"

manifest_path(minor::VersionNumber) = joinpath(BIN, "Manifest-v$(series(minor)).toml")

# The Julia minors bin/ supports: released versions only (a prerelease bundles
# what it will bundle, not what it does), matching the project's own bound.
function supported_minors()
    project = Pkg.TOML.parsefile(joinpath(BIN, "Project.toml"))
    # Project.toml compat is semver syntax ("^1.9"), not registry syntax
    bound = Pkg.Types.semver_spec(get(get(project, "compat", Dict()), "julia", "*"))
    minors = Set{VersionNumber}()
    for v in Registries.JULIA_VERSIONS
        isempty(v.prerelease) && v in bound || continue
        push!(minors, VersionNumber(v.major, v.minor))
    end
    return sort!(collect(minors))
end

# Is `julia +X.Y` a channel juliaup can run? Probing beats parsing `juliaup
# status`, and it is the same question the update itself will ask.
function channel_available(minor::VersionNumber)
    channel = "+$(series(minor))"
    success(pipeline(`julia $channel --version`; stdout = devnull, stderr = devnull))
end

function update_manifest(minor::VersionNumber)
    channel = "+$(series(minor))"
    # the child writes straight to fd 1, so anything still buffered here would
    # otherwise land after output from a subprocess that ran later
    flush(stdout)
    run(`julia $channel --project=$BIN -e "using Pkg; Pkg.update()"`)
    flush(stdout)
end

# Which supported minors have no manifest of their own? This is the staleness
# that no amount of testing catches, because the symptom is a Julia quietly
# resolving against another Julia's pins rather than anything failing.
function report_missing(minors::Vector{VersionNumber})
    absent = filter(m -> !isfile(manifest_path(m)), minors)
    if isempty(absent)
        println("every supported Julia minor has a manifest")
        return true
    end
    println("no manifest for: ", join(series.(absent), ", "))
    println("\nRun bin/update_manifests.jl on a machine with those Julias installed:\n",
            "    juliaup add ", join(series.(absent), " "), "\n",
            "    julia bin/update_manifests.jl")
    # The snapshot matrix in the workflow is spelled out rather than derived, so
    # a new minor needs adding there too -- and nothing else would say so: the
    # matrix would keep passing on the versions it does list.
    println("\nThen add ", join(("'" * series(m) * "'" for m in absent), ", "),
            " to the `snapshot` matrix in .github/workflows/stdlibs.yml.")
    return false
end

function main(args::Vector{String} = ARGS)
    check_only = "--check" in args
    isempty(setdiff(args, ["--check"])) ||
        error("usage: update_manifests.jl [--check]")
    minors = supported_minors()
    isempty(minors) && error("no released Julia matches bin/Project.toml's julia compat")
    println("supported Julia minors: ", join(series.(minors), ", "), "\n")
    check_only && return report_missing(minors)

    missing_channels, changed, created = VersionNumber[], VersionNumber[], VersionNumber[]
    for minor in minors
        path = manifest_path(minor)
        if !channel_available(minor)
            push!(missing_channels, minor)
            println("julia $(series(minor)): channel not installed, skipping")
            continue
        end
        # A Julia with no manifest of its own falls back to Manifest.toml, which
        # is the 1.9 symlink -- so it would overwrite 1.9's pins with its own.
        # Giving it an empty manifest first makes it resolve into the right file.
        is_new = !isfile(path)
        is_new && touch(path)
        before = read(path, String)
        println("\n=== julia $(series(minor)) ===")
        try
            update_manifest(minor)
        catch
            # an empty file left behind would be committed as a manifest and
            # then quietly used as one, which is worse than having none
            is_new && rm(path; force = true)
            rethrow()
        end
        if is_new
            push!(created, minor)
        elseif read(path, String) != before
            push!(changed, minor)
        end
    end

    println()
    # The 1.9 manifest is reached through a symlink, since 1.9 does not
    # understand version-specific manifest files. Pkg writing through it should
    # leave the link alone, but a manifest silently turned into a regular file
    # would strand 1.9 on stale pins, so say so rather than assume.
    link = joinpath(BIN, "Manifest.toml")
    ispath(link) && !islink(link) &&
        println("WARNING: $link is no longer a symlink -- restore it with\n",
                "    ln -sf Manifest-v1.9.toml $link")
    isempty(created) || println("created: ", join(("Manifest-v$(series(m)).toml" for m in created), ", "))
    if isempty(changed)
        println("all manifests already current")
    else
        println("updated: ", join(("Manifest-v$(series(m)).toml" for m in changed), ", "))
        println("\nReview and commit the changes, and run the bin/ tests:\n",
                "    julia --project=bin bin/test/runtests.jl")
    end
    isempty(missing_channels) ||
        println("\nNot checked -- these channels are not installed:\n",
                "    juliaup add ", join(series.(missing_channels), " "))
    return isempty(missing_channels)
end

main() || exit(1)
