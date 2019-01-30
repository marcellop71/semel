{
  description = "Semel – simplicial methods library";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs = { self, nixpkgs }:
    let
      supportedSystems = [ "x86_64-linux" "aarch64-linux" ];
      forAllSystems = nixpkgs.lib.genAttrs supportedSystems;
    in
    {
      packages = forAllSystems (system:
        let
          pkgs = nixpkgs.legacyPackages.${system};

          machine = pkgs.stdenv.hostPlatform.uname.processor;

          deps = [
            pkgs.zlog
            pkgs.judy
            pkgs.flint3
            pkgs.gsl
            pkgs.gmp
            pkgs.mpfr
            pkgs.qhull
          ];
        in
        {
          default = pkgs.stdenv.mkDerivation {
            pname = "semel";
            version = "0.1.0";

            src = ./.;

            nativeBuildInputs = [ pkgs.cmake ];
            buildInputs = deps;

            cmakeFlags = [
              "-DB_PROJECT_NAME=semel"
              "-DB_VERSION_MAJOR=0"
              "-DB_VERSION_MINOR=1"
              "-DB_VERSION_PATCH=0"
              "-DB_MACHINE=${machine}"
              "-DZLOG_CONF_PATH=${placeholder "out"}/etc/semel/zlog.conf"
            ];

            postInstall = ''
              mkdir -p $out/etc/semel
              cp ${./config/zlog.conf} $out/etc/semel/zlog.conf

              for h in ${./include}/*.h; do
                install -Dm644 "$h" "$out/include/semel/$(basename "$h")"
              done
            '';
          };
        }
      );

      devShells = forAllSystems (system:
        let
          pkgs = nixpkgs.legacyPackages.${system};
        in
        {
          default = pkgs.mkShell {
            inputsFrom = [ self.packages.${system}.default ];

            packages = [
              pkgs.python3
              pkgs.python3Packages.numpy
              pkgs.python3Packages.matplotlib
            ];
          };
        }
      );
    };
}
