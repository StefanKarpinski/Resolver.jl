using Pkg, UUIDs
@show deps = Pkg.dependencies()
@show pkg = deps[UUID("34da2185-b29b-5c13-b0c7-acf172513d20")]
@show compat_version = pkg.version
# exit(compat_version == v"4.0.0" ? 1 : 0)
