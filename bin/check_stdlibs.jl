#!/usr/bin/env julia
#
# Check the pinned stdlib snapshot against the Julia running this script.
#
#     julia --project=bin bin/check_stdlibs.jl
#
# `bin/` pins HistoricalStdlibVersions in its manifests, and that package is a
# snapshot: it records what each Julia release bundles. A pin left behind by new
# Julia releases does not fail loudly, it answers wrongly -- a release newer than
# the pin is not missing from the data, it is described by an older stanza. So
# the resolver reports confidently about a Julia that never existed, and the
# report looks exactly like a real conflict. That is issue #103: the pin sat at
# 2.0.2, which predates 1.12, so 1.12.x was credited with 1.11's OpenSSL_jll
# 3.0.15 rather than the 3.5.x it ships, and OpenSSL_CLI_jll would not resolve
# on it at all.
#
# The running Julia carries the answer in its own stdlib directory, so that is
# what the snapshot is checked against. Note what this does *not* do: it never
# asks whether a newer HistoricalStdlibVersions exists, so it cannot fail merely
# because time passed -- only because the pin actually misdescribes a Julia that
# is present. Keeping the pin current is bin/update_manifests.jl's job.
#
# Coverage is one Julia per run, so .github/workflows/stdlibs.yml runs it across
# the range the bin/ manifests support. Prereleases are skipped: see below.
# Exits nonzero on drift.

push!(empty!(LOAD_PATH), @__DIR__)

include("Registries.jl")

function check_stdlibs(io::IO = stdout)
    # Only a released Julia has a settled answer to "what do you bundle". A
    # development build is a moving target -- `1.14.0-DEV.2617` is master at one
    # commit, while the snapshot's 1.14.0 stanza is master at another -- so the
    # two disagree for no reason anybody could fix, and on this machine's
    # nightly they disagree about eighteen stdlibs. Released versions are also
    # the only ones the bin/ manifests target, so skipping the rest costs the
    # check nothing it was meant to catch.
    if !isempty(VERSION.prerelease)
        println(io, "julia $VERSION is a prerelease: skipping the stdlib snapshot check")
        return true
    end
    bundled = host_stdlibs()
    snapshot = stdlib_snapshot(VERSION)
    names = Dict(uuid => info.name for (uuid, info) in snapshot)
    # The upgradable stdlibs are absent from the snapshot by design -- Julia
    # bundles them without pinning them, so registry versions compete with the
    # bundled one -- and their absence is not staleness.
    pinned = Dict(uuid => info.version
                  for (uuid, info) in snapshot
                  if info.version !== nothing &&
                     uuid ∉ UPGRADABLE_STDLIBS_UUIDS)
    absent = sort!([names[uuid] for uuid in setdiff(keys(pinned), keys(bundled))])
    drifted = sort!([(names[uuid], version, bundled[uuid])
                     for (uuid, version) in pinned
                     if haskey(bundled, uuid) && bundled[uuid] != version])
    println(io, "julia $VERSION: $(length(pinned)) stdlibs pinned by the snapshot, ",
                "$(length(bundled)) bundled")
    if isempty(absent) && isempty(drifted)
        println(io, "the pinned stdlib snapshot describes this Julia")
        return true
    end
    for name in absent
        println(io, "  missing: $name is pinned by the snapshot but not bundled here")
    end
    for (name, claimed, actual) in drifted
        println(io, "  drifted: $name -- snapshot says $claimed, this Julia ships $actual")
    end
    println(io, """

        The stdlib snapshot does not describe julia $VERSION. Refresh the pin:

            julia --project=bin bin/update_manifests.jl

        and commit the bin/Manifest-*.toml changes.""")
    return false
end

check_stdlibs() || exit(1)
