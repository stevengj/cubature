/* Copyright (c) 2007-2010 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
 */

/* simple implementation of red-black trees optimized for use with
   with scubature.c, based on code I originally wrote for my NLopt
   library (http://ab-initio.mit.edu/nlopt).  -- Steven G. Johnson */

/* USAGE: 

   I've hacked my original redblack implementation to make it
   usable only from one file, without polluting the global namespace,
   with a hardcoded key type and comparison function.  To use
   it more generally, you probably want the redblack.c and redblack.h
   files from my NLopt library (see above).   The original NLopt version
   is also more full-featured (more search functions, deletion, etc.).

   Here:

   #include <stdlib.h>

   typedef ....whatever.... rb_key;
   typedef ....whatever.... rb_key_data;
   static void rb_destroy_key(rb_key k) { ... }
   static int rb_compare(const rb_key k1, const rb_key k2,
                         rb_key_data data) { ... return +1/0/-1 ... }

   #include "redblack.c"

   ....

   Functions ending with an underscore (_) are internal functions not
   designed to be called outside of redblack.c
*/

typedef enum { RED, BLACK } rb_color;
typedef struct rb_node_s {
     struct rb_node_s *p, *r, *l; /* parent, right, left */
     rb_key k; /* key (and data) */
     rb_color c;
} rb_node;

typedef struct {
     rb_node *root;
     rb_key_data data; /* shared data for comparison function */
     int N; /* number of nodes */
} rb_tree;

/* it is convenient to use an explicit node for NULL nodes ... we need
   to be careful never to change this node indirectly via one of our
   pointers!  */
static rb_node rb_nil = {&rb_nil, &rb_nil, &rb_nil, 0, BLACK};
#define NIL (&rb_nil)

static void rb_tree_init(rb_tree *t, rb_key_data data) {
     t->data = data;
     t->root = NIL;
     t->N = 0;
}

static void rb_destroy_(rb_node *n)
{
     if (n != NIL) {
	  rb_destroy_(n->l); rb_destroy_(n->r);
	  rb_destroy_key(n->k);
	  free(n);
     }
}

static void rb_tree_destroy(rb_tree *t)
{
     rb_destroy_(t->root);
     t->root = NIL;
}

static void rb_rotate_left_(rb_node *p, rb_tree *t)
{
     rb_node *n = p->r; /* must be non-NIL */
     p->r = n->l;
     n->l = p;
     if (p->p != NIL) {
	  if (p == p->p->l) p->p->l = n;
	  else p->p->r = n;
     }
     else
	  t->root = n;
     n->p = p->p;
     p->p = n;
     if (p->r != NIL) p->r->p = p;
}

static void rb_rotate_right_(rb_node *p, rb_tree *t)
{
     rb_node *n = p->l; /* must be non-NIL */
     p->l = n->r;
     n->r = p;
     if (p->p != NIL) {
	  if (p == p->p->l) p->p->l = n;
	  else p->p->r = n;
     }
     else
	  t->root = n;
     n->p = p->p;
     p->p = n;
     if (p->l != NIL) p->l->p = p;
}

static void rb_insert_node_(rb_tree *t, rb_node *n)
{
     rb_key_data data = t->data;
     rb_key k = n->k;
     rb_node *p = t->root;
     n->c = RED;
     n->p = n->l = n->r = NIL;
     t->N++;
     if (p == NIL) {
	  t->root = n;
	  n->c = BLACK;
	  return;
     }
     /* insert (RED) node into tree */
     while (1) {
	  if (rb_compare(k, p->k, data) <= 0) { /* k <= p->k */
	       if (p->l != NIL)
		    p = p->l;
	       else {
		    p->l = n;
		    n->p = p;
		    break;
	       }
	  }
	  else {
	       if (p->r != NIL)
		    p = p->r;
	       else {
		    p->r = n;
		    n->p = p;
		    break;
	       }
	  }
     }
 fixtree:
     if (n->p->c == RED) { /* red cannot have red child */
	  rb_node *u = p == p->p->l ? p->p->r : p->p->l;
	  if (u != NIL && u->c == RED) {
	       p->c = u->c = BLACK;
	       n = p->p;
	       if ((p = n->p) != NIL) {
		    n->c = RED;
		    goto fixtree;
	       }
	  }
	  else {
	       if (n == p->r && p == p->p->l) {
		    rb_rotate_left_(p, t);
		    p = n; n = n->l;
	       }
	       else if (n == p->l && p == p->p->r) {
		    rb_rotate_right_(p, t);
		    p = n; n = n->r;
	       }
	       p->c = BLACK;
	       p->p->c = RED;
	       if (n == p->l && p == p->p->l)
		    rb_rotate_right_(p->p, t);
	       else if (n == p->r && p == p->p->r)
		    rb_rotate_left_(p->p, t);
	  }
	      
     }
}

static rb_node *rb_tree_insert(rb_tree *t, rb_key k)
{
     rb_node *n = (rb_node *) malloc(sizeof(rb_node));
     if (!n) return NULL;
     n->k = k;
     rb_insert_node_(t, n);
     return n;
}

static rb_node *rb_tree_find(rb_tree *t, rb_key k)
{
     rb_key_data data = t->data;
     rb_node *p = t->root;
     while (p != NIL) {
	  int comp = rb_compare(k, p->k, data);
	  if (!comp) return p;
	  p = comp <= 0 ? p->l : p->r;
     }
     return NULL;
}
