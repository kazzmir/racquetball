import os

env = Environment(ENV = os.environ)

env.VariantDir('build', 'src')
env.ParseConfig('pkg-config allegro-5 --cflags --libs')
env.Append(CCFLAGS = ['-g3'])
source = Split("""build/main.cpp""")

env.Program('racquetball', source)
