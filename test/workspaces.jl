# Integration tests for Julia 1.12+ workspace support in bin/resolve.jl.
#
# Runnable standalone from the repo root (needs the General registry
# installed and takes a few minutes — each case shells out to a full
# registry resolve):
#
#     julia +release --project=. test/workspaces.jl
#
# Julia 1.12+ resolves all projects in a workspace into a single shared
# manifest. These tests need >= 1.12 for the feature, and -- since we generate
# a full manifest -- a released Julia, so the resolution target matches the
# host (a prerelease host has no matching registered release to target).

using Test
using Resolver
import TOML

# Run bin/resolve.jl on a workspace; return (success, stdout, stderr strings).
function resolve_workspace_manifest(ws::AbstractString; flags::Cmd = `--print-manifest`)
    julia = Base.julia_cmd()[1]
    project = pkgdir(Resolver, "bin")
    script = joinpath(project, "resolve.jl")
    julia_version = string(Int(VERSION.major), ".", VERSION.minor)
    run(`$julia --project=$project -e 'import Pkg; Pkg.instantiate()'`)
    out = IOBuffer()
    err = IOBuffer()
    ok = success(pipeline(
        `$julia --project=$project $script $ws $flags --julia=$julia_version`;
        stdout = out, stderr = err,
    ))
    return ok, String(take!(out)), String(take!(err))
end

@testset "workspaces" begin
    if VERSION < v"1.12" || !isempty(VERSION.prerelease)
        return
    end

    @testset "interdependent workspace packages" begin
        # a workspace where the root and two members are all packages:
        #  - the root depends on a member (the monorepo shape)
        #  - the second member depends on the root and on the first member
        #  - the root and first member each have a registry dep of their own
        ws = mktempdir()
        mkdir(joinpath(ws, "member"))
        mkdir(joinpath(ws, "member2"))
        write(joinpath(ws, "Project.toml"), """
            name = "WorkspaceRoot"
            uuid = "11111111-1111-1111-1111-111111111111"
            version = "0.1.0"

            [deps]
            Example = "7876af07-990d-54b4-ab0e-23690620f79a"
            WorkspaceMember = "22222222-2222-2222-2222-222222222222"

            [workspace]
            projects = ["member", "member2"]
            """)
        write(joinpath(ws, "member", "Project.toml"), """
            name = "WorkspaceMember"
            uuid = "22222222-2222-2222-2222-222222222222"
            version = "0.2.0"

            [deps]
            JSON = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
            """)
        write(joinpath(ws, "member2", "Project.toml"), """
            name = "WorkspaceExtra"
            uuid = "33333333-3333-3333-3333-333333333333"
            version = "0.3.0"

            [deps]
            WorkspaceMember = "22222222-2222-2222-2222-222222222222"
            WorkspaceRoot = "11111111-1111-1111-1111-111111111111"
            """)

        ok, manifest, err = resolve_workspace_manifest(ws)
        ok || print(stderr, err)
        @test ok
        m = TOML.parse(manifest)

        # all three workspace packages are path entries relative to the root,
        # with their local versions and inter-workspace deps preserved
        root = only(m["deps"]["WorkspaceRoot"])
        @test root["path"] == "."
        @test root["version"] == "0.1.0"
        @test sort(root["deps"]) == ["Example", "WorkspaceMember"]

        member = only(m["deps"]["WorkspaceMember"])
        @test member["path"] == "member"
        @test member["version"] == "0.2.0"
        @test member["deps"] == ["JSON"]

        extra = only(m["deps"]["WorkspaceExtra"])
        @test extra["path"] == "member2"
        @test extra["version"] == "0.3.0"
        @test sort(extra["deps"]) == ["WorkspaceMember", "WorkspaceRoot"]

        # the workspace packages' registry deps resolve as registry entries
        for name in ("Example", "JSON")
            entry = only(m["deps"][name])
            @test haskey(entry, "git-tree-sha1")
            @test !haskey(entry, "path")
        end
    end

    @testset "member shadowing a registered package" begin
        # a workspace member that is a local development copy of a registered
        # package (JSON), at a version that does not exist in the registry,
        # with the root's compat targeting that local version. Pkg treats
        # workspace packages as fixed path entries and never consults the
        # registry for them; check that we do the same.
        ws = mktempdir()
        mkdir(joinpath(ws, "json"))
        write(joinpath(ws, "Project.toml"), """
            name = "WorkspaceRoot"
            uuid = "11111111-1111-1111-1111-111111111111"
            version = "0.1.0"

            [deps]
            JSON = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"

            [compat]
            JSON = "99"

            [workspace]
            projects = ["json"]
            """)
        write(joinpath(ws, "json", "Project.toml"), """
            name = "JSON"
            uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
            version = "99.0.0"

            [deps]
            Example = "7876af07-990d-54b4-ab0e-23690620f79a"
            """)

        ok, manifest, err = resolve_workspace_manifest(ws)
        ok || print(stderr, err)
        @test ok
        m = TOML.parse(manifest)

        # the local JSON is a path entry fixed at its local version, and the
        # registry version of JSON (and its dependency tree: Parsers etc.)
        # appears nowhere in the manifest
        json = only(m["deps"]["JSON"])
        @test json["path"] == "json"
        @test json["version"] == "99.0.0"
        @test json["deps"] == ["Example"]
        @test sort(collect(keys(m["deps"]))) == ["Example", "JSON", "WorkspaceRoot"]
    end

    @testset "bare workspace projects" begin
        # (a) a bare sub-environment member (test/docs style: deps but no
        # name/uuid): its deps are folded into resolution, but it is not
        # itself a manifest entry -- only the root gets a path entry
        ws = mktempdir()
        mkdir(joinpath(ws, "test"))
        write(joinpath(ws, "Project.toml"), """
            name = "WorkspaceRoot"
            uuid = "11111111-1111-1111-1111-111111111111"
            version = "0.1.0"

            [deps]
            Example = "7876af07-990d-54b4-ab0e-23690620f79a"

            [workspace]
            projects = ["test"]
            """)
        write(joinpath(ws, "test", "Project.toml"), """
            [deps]
            JSON = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
            WorkspaceRoot = "11111111-1111-1111-1111-111111111111"
            """)

        ok, manifest, err = resolve_workspace_manifest(ws)
        ok || print(stderr, err)
        @test ok
        m = TOML.parse(manifest)
        root = only(m["deps"]["WorkspaceRoot"])
        @test root["path"] == "."
        @test haskey(only(m["deps"]["JSON"]), "git-tree-sha1")
        @test haskey(m["deps"], "Example")
        # exactly one path entry: the root (the bare member has none)
        @test sum(haskey(only(e), "path") for e in values(m["deps"])) == 1

        # (b) a bare root (no name/uuid) whose member is a package: only the
        # member gets a path entry
        ws = mktempdir()
        mkdir(joinpath(ws, "member"))
        write(joinpath(ws, "Project.toml"), """
            [deps]
            Example = "7876af07-990d-54b4-ab0e-23690620f79a"

            [workspace]
            projects = ["member"]
            """)
        write(joinpath(ws, "member", "Project.toml"), """
            name = "WorkspaceMember"
            uuid = "22222222-2222-2222-2222-222222222222"
            version = "0.2.0"

            [deps]
            JSON = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
            """)

        ok, manifest, err = resolve_workspace_manifest(ws)
        ok || print(stderr, err)
        @test ok
        m = TOML.parse(manifest)
        member = only(m["deps"]["WorkspaceMember"])
        @test member["path"] == "member"
        @test member["deps"] == ["JSON"]
        @test haskey(m["deps"], "Example")
        @test haskey(m["deps"], "JSON")
        @test sum(haskey(only(e), "path") for e in values(m["deps"])) == 1
    end

    @testset "registry package depending on a workspace package" begin
        # Conda (registry) depends on JSON; the workspace shadows JSON with a
        # local dev copy. The local copy is fixed: resolution must use its
        # local version and deps -- never the registry's JSON -- and
        # dependents' compat must be enforced against the fixed version,
        # exactly as Pkg does. Conda is used because its compat on JSON is a
        # bounded range: it admits 1.x but can never admit 99.x.
        function make_shadow_ws(json_version)
            ws = mktempdir()
            mkdir(joinpath(ws, "json"))
            write(joinpath(ws, "Project.toml"), """
                name = "WorkspaceRoot"
                uuid = "11111111-1111-1111-1111-111111111111"
                version = "0.1.0"

                [deps]
                Conda = "8f4d0f93-b110-5947-807f-2305c1781a2d"

                [workspace]
                projects = ["json"]
                """)
            write(joinpath(ws, "json", "Project.toml"), """
                name = "JSON"
                uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
                version = "$json_version"

                [deps]
                Example = "7876af07-990d-54b4-ab0e-23690620f79a"
                """)
            return ws
        end

        # (a) local version admitted by Conda's compat: resolves against the
        # fixed local copy; the registry JSON's dependency tree stays out
        ws = make_shadow_ws("1.99.99")
        ok, out, err = resolve_workspace_manifest(ws; flags = `--print-versions`)
        ok || print(stderr, err)
        @test ok
        @test occursin(r"(?m)^682c06a0-de6a-54ab-a142-c8b1cf79cde6 JSON\s+1\.99\.99$", out)
        @test !occursin("Parsers", out)

        ok, manifest, err = resolve_workspace_manifest(ws)
        ok || print(stderr, err)
        @test ok
        m = TOML.parse(manifest)
        json = only(m["deps"]["JSON"])
        @test json["path"] == "json"
        @test json["version"] == "1.99.99"
        @test json["deps"] == ["Example"]
        conda = only(m["deps"]["Conda"])
        @test haskey(conda, "git-tree-sha1")
        @test "JSON" in conda["deps"]
        @test haskey(m["deps"], "Example")
        @test !haskey(m["deps"], "Parsers")

        # (b) local version excluded by Conda's compat: unsatisfiable, as Pkg
        ws = make_shadow_ws("99.0.0")
        ok, out, err = resolve_workspace_manifest(ws)
        @test !ok
        @test occursin("Unresolvable", err)   # resolve.jl now prints a diagnosis
    end

    @testset "compat excluding a workspace package's local version" begin
        # compat naming a workspace package must admit its fixed local
        # version, matching Pkg's enforcement for fixed packages
        ws = mktempdir()
        mkdir(joinpath(ws, "json"))
        write(joinpath(ws, "Project.toml"), """
            name = "WorkspaceRoot"
            uuid = "11111111-1111-1111-1111-111111111111"
            version = "0.1.0"

            [deps]
            JSON = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"

            [compat]
            JSON = "1"

            [workspace]
            projects = ["json"]
            """)
        write(joinpath(ws, "json", "Project.toml"), """
            name = "JSON"
            uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
            version = "99.0.0"
            """)

        ok, out, err = resolve_workspace_manifest(ws)
        @test !ok
        @test occursin("excludes its local version", err)
    end
end
