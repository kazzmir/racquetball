import os

env = Environment(ENV = os.environ)

env.VariantDir('build', 'src')
source = Split("""build/main.cpp""")

env.Program('racquetball', source)
