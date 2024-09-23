from triqs_dft_tools.converters import Wannier90Converter

Converter = Wannier90Converter(seedname="nsp", w90zero=3e-6)
Converter.convert_dft_input()
