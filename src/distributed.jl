function init_workers(nprocs=35)
    rmprocs(workers());
    addprocs(nprocs)
    exp = quote
        @everywhere begin
            import Pkg
        Pkg.activate(".")
            using Lensman
        end
    end
    @suppress_err eval(macroexpand(Distributed, exp))
end
