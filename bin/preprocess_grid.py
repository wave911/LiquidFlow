def subfinder(items, pattern):
    res = list(filter(lambda x: x in pattern, items))
    if (len(res) == len(pattern)):
    	return res
    else:
    	return None

#flags = tf.app.flags
#FLAGS = flags.FLAGS

#flags.DEFINE_string('mesh_file', '', "mesh file (.data)")

f = open('../mesh/Mesh_box.dat','a+')
f.seek(0,0)
pntCount, elemTotalCount = f.readline().split(' ')

points = list()
triangles = list()
elements = dict()
triangle_type = ''
idx = 0
pts = 0
for line in f:
	if (pts < int(pntCount)):
		pts += 1
		continue
	num, type, pnts = line.split(' ', 2)
	if (type.startswith('20')):
		triangle_type = type
		pts_list = pnts.split(' ')[:-1]
		triangles.append(list(map(int, pts_list)))
	if (type.startswith('30')):
		pts_list = pnts.split(' ')[:-1]
		elements[idx] = list(map(int, pts_list))
		idx += 1


for t in triangles:
	for k, v in elements.items():
		res = subfinder(sorted(v), sorted(t))
		if (res is not None):
			f.write(str(k) + ' ' + '2' + triangle_type + ' ' + str(t).strip('[]') + '\n')
f.close()
