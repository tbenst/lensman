{ pkgs ? import <nixpkgs> {} }:

let 
  julia-env = pkgs.callPackage /Computer/packages/julia-env.nix {};
in
pkgs.mkShell {
  buildInputs = with pkgs; [
    conda
    vscode
    julia-env
  ];

  LD_LIBRARY_PATH="/run/opengl-driver/lib:/run/opengl-driver-32/lib";
}
