
def mergeSortInversions(arr):
    if len(arr) == 1:
        return arr, 0
    else:
        a = arr[:len(arr)//2]
        b = arr[len(arr)//2:]
        a, ai = mergeSortInversions(a)
        b, bi = mergeSortInversions(b)
        c = []
        i = 0
        j = 0
        inversions = 0 + ai + bi
    while i < len(a) and j < len(b):
        if a[i] <= b[j]:
            c.append(a[i])
            i += 1
        else:
            c.append(b[j])
            j += 1
            inversions += (len(a)-i)
    c += a[i:]
    c += b[j:]
    return c, inversions

def Inversions(arr):
    b,c=mergeSortInversions(arr)
    return c


def BraidOrbit(word,rels):
 def pattern_match (L, i, X, l):
  for ind in range(l):
   if L[i+ind] != X[ind]:
    return False
  return True

 l=len(word)
 words = set(tuple(word))
 test_words = [ tuple(word) ]

 rels = rels + [ [b,a] for a,b in rels ]
 rels = [ [tuple(a), tuple(b), len(a) ]  for a,b in rels ]

 loop_ind = 0
 list_len = 1
 yield word
 while loop_ind < list_len:
  test_word = test_words[loop_ind]
  loop_ind += 1
  for rel in rels:
   left = rel[0]
   right = rel[1]
   rel_l = rel[2]
   for i in range(l-rel_l+1):
    if pattern_match(test_word, i, left, rel_l):
     new_word=test_word[:i]+right+test_word[i+rel_l:]
     if new_word not in words:
      words.add(new_word)
      test_words.append(new_word)
      list_len+=1
      yield new_word
 return

def coxeter_braid_orbit(coxeter_group, word):
 word=list(word)
 braid_rels=coxeter_group.braid_relations()
 return BraidOrbit(word, braid_rels)


def braid_orbit_generator(wg,word,num,inversions):
 gen=coxeter_braid_orbit(wg,word)
 n=0
 for el in gen:
  y=True
  if inversions>=0 and Inversions(el)!=inversions:
   y=False
  if num>0 and n>=num:
   break
  if y:
   n+=1
   yield el
 return

