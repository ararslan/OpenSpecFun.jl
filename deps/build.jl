using BinDeps

@BinDeps.setup

osf = library_dependency("openspecfun",
                         aliases=["libopenspecfun", "libopenspecfun.1", "libopenspecfun.1.3"])

# if is_apple()
#     using Homebrew
#     provides(Homebrew.HB, "ararslan/pints/openspecfun", osf, os=:Darwin)
# else
    if isdir(srcdir(osf))
        rm(srcdir(osf), recursive=true)
        mkdir(srcdir(osf))
    end

    if isdir(BinDeps.downloadsdir(osf))
        rm(BinDeps.downloadsdir(osf), recursive=true)
        mkdir(BinDeps.downloadsdir(osf))
    end

    vers = v"0.5.3"

    provides(Sources, URI("https://github.com/JuliaLang/openspecfun/archive/v$vers.tar.gz"), osf,
            unpacked_dir="")

    provides(BuildProcess, (@build_steps begin
        GetSources(osf)
        CreateDirectory(joinpath(BinDeps.builddir(osf), "openspecfun"))
        @build_steps begin
            ChangeDirectory(joinpath(BinDeps.builddir(osf), "openspecfun"))
            FileRule(joinpath(libdir(osf), "libopenspecfun." * BinDeps.shlib_ext), @build_steps begin
                CreateDirectory(libdir(osf))
                `pwd`
                `ls -R`
                `make install`
            end)
        end
    end), osf)
# end

@BinDeps.install Dict(:openspecfun => :openspecfun)
