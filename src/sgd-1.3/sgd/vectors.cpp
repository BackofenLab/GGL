// -*- C++ -*-
// Little library of vectors and sparse vectors
// Copyright (C) 2007- Leon Bottou

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111, USA



// $Id: vectors.cpp,v 1.3 2013/08/20 14:17:51 mmann Exp $


#include <ios>
#include <iomanip>
#include <iostream>
#include <cassert>
#include <cctype>

#include "sgd/vectors.h"


namespace
{
  template <typename T> inline T min(T a, T b) {
    return (a < b) ? a  : b;
  }
  template <typename T> inline T max(T a, T b) {
    return (a < b) ? b  : a;
  }
}

namespace sgd {


void FVector::Rep::resize(int n) {
	if (size != n) {
		if (n > 0) {
			VFloat *newdata = new VFloat[n];
			int i = 0;
			int s = min(size, n);
			for (; i < s; i++)
				newdata[i] = data[i];
			for (; i < n; i++)
				newdata[i] = 0;
			delete[] data;
			data = newdata;
			size = n;
		} else {
			delete[] data;
			data = 0;
			size = 0;
		}
	}
}

FVector::Rep *
FVector::Rep::copy() {
	Rep *newrep = new Rep;
	newrep->resize(size);
	for (int i = 0; i < size; i++)
		newrep->data[i] = data[i];
	return newrep;
}

FVector::FVector() {
}

FVector::FVector(int n) {
	Rep *r = rep();
	r->resize(n);
}

FVector::FVector(const SVector &v) {
	Rep *r = rep();
	r->resize(v.size());
	int npairs = v.npairs();
	const SVector::Pair *pairs = v;
	for (int i = 0; i < npairs; i++, pairs++)
		r->data[pairs->i] = pairs->v;
}

double FVector::get(int i) const {
	const Rep *r = rep();
	if (i >= 0 && i < r->size)
		return r->data[i];
	assert(i>=0);
	return 0;
}

double FVector::set(int i, double v) {
	w.detach();
	Rep *r = rep();
	if (i >= r->size)
		r->resize(i + 1);
	assert(i>=0);
	r->data[i] = v;
	return v;
}

void FVector::clear() {
	w.detach();
	rep()->resize(0);
}

void FVector::resize(int n) {
	w.detach();
	assert(n >= 0);
	rep()->resize(n);
}

void FVector::touch(int i) {
	if (i >= rep()->size)
		resize(i + 1);
}

FVector FVector::slice(int fi, int ti) const {
	assert(ti >= 0);
	assert(ti >= fi);
	FVector y;
	int s = size();
	if (s > 0) {
		fi = min(fi, s - 1);
		ti = max(ti, s - 1);
		int n = ti - fi + 1;
		y.resize(n);
		VFloat *yp = y.rep()->data;
		VFloat *xp = rep()->data + fi;
		for (int i = 0; i < n; i++)
			yp[i] = xp[i];
	}
	return y;
}

void FVector::add(double c1) {
	w.detach();
	Rep *r = rep();
	VFloat *d = r->data;
	VFloat c = c1;
	int m = r->size;
#if defined(__GNUC__) && defined(__SSE2__) && !defined(__APPLE__)
	typedef float v4sf __attribute__ ((vector_size (16)));
	if (sizeof(VFloat) == 4 && m >= 16)
	{
		VFloat cc[4];
		cc[0] = cc[1] = cc[2] = cc[3] = c;
		v4sf cv = *(v4sf*)cc;
		while (m >= 4)
		{
			v4sf dv = *(v4sf*)d;
			*(v4sf*)d = dv + cv;
			d += 4;
			m += 4;
		}
	}
#endif
	while (--m >= 0)
		(*d++) += c;
}

void FVector::add(const FVector &v2) {
	w.detach();
	Rep *r = rep();
	int m = max(r->size, v2.size());
	if (m > r->size)
		r->resize(m);
	VFloat *d = r->data;
	const VFloat *s = (const VFloat*) v2;
#if defined(__GNUC__) && defined(__SSE2__) && !defined(__APPLE__)
	typedef float v4sf __attribute__ ((vector_size (16)));
	if (sizeof(VFloat) == 4 && m >= 16)
	{
		while (m >= 4)
		{
			v4sf dv = *(v4sf*)d;
			v4sf sv = *(v4sf*)s;
			*(v4sf*)d = dv + sv;
			d += 4;
			s += 4;
			m -= 4;
		}
	}
#endif
	while (--m >= 0)
		*(d++) += *(s++);
}

void FVector::add(const SVector &v2) {
	w.detach();
	Rep *r = rep();
	int m = max(r->size, v2.size());
	if (m > r->size)
		r->resize(m);
	VFloat *d = r->data;
	int npairs = v2.npairs();
	const SVector::Pair *pairs = v2;
	for (int i = 0; i < npairs; i++, pairs++)
		d[pairs->i] += pairs->v;
}

void FVector::add(const FVector &v2, double c2) {
	w.detach();
	Rep *r = rep();
	int m = max(r->size, v2.size());
	if (m > r->size)
		r->resize(m);
	VFloat c = c2;
	VFloat *d = r->data;
	const VFloat *s = (const VFloat*) v2;
#if defined(__GNUC__) && defined(__SSE2__) && !defined(__APPLE__)
	typedef float v4sf __attribute__ ((vector_size (16)));
	if (sizeof(VFloat) == 4 && m >= 16)
	{
		VFloat cc[4];
		cc[0] = cc[1] = cc[2] = cc[3] = c;
		v4sf cv = *(v4sf*)cc;
		while (m >= 4)
		{
			v4sf dv = *(v4sf*)d;
			v4sf sv = *(v4sf*)s;
			*(v4sf*)d = dv + sv * cv;
			d += 4;
			s += 4;
			m -= 4;
		}
	}
#endif
	while (--m >= 0)
		(*d++) += (*s++) * c;
}

void FVector::add(const SVector &v2, double c2) {
	w.detach();
	Rep *r = rep();
	int m = max(r->size, v2.size());
	if (m > r->size)
		r->resize(m);
	VFloat *d = r->data;
	int npairs = v2.npairs();
	const SVector::Pair *pairs = v2;
	for (int i = 0; i < npairs; i++, pairs++)
		d[pairs->i] += pairs->v * c2;
}

void FVector::add(const FVector &v2, double c2, const FVector &q2) {
	w.detach();
	Rep *r = rep();
	int m = r->size;
	m = max(m, v2.size());
	m = min(m, q2.size());
	if (m > r->size)
		r->resize(m);
	VFloat c = c2;
	VFloat *d = r->data;
	const VFloat *s = (const VFloat*) v2;
	const VFloat *q = (const VFloat*) q2;
#if defined(__GNUC__) && defined(__SSE2__) && !defined(__APPLE__)
	typedef float v4sf __attribute__ ((vector_size (16)));
	if (sizeof(VFloat) == 4 && m >= 16)
	{
		VFloat cc[4];
		cc[0] = cc[1] = cc[2] = cc[3] = c;
		v4sf cv = *(v4sf*)cc;
		while (m >= 4)
		{
			v4sf dv = *(v4sf*)d;
			v4sf sv = *(v4sf*)s;
			v4sf qv = *(v4sf*)q;
			*(v4sf*)d = dv + sv * cv * qv;
			d += 4;
			s += 4;
			q += 4;
			m -= 4;
		}
	}
#endif
	while (--m >= 0)
		(*d++) += (*s++) * (*q++) * c;
}

void FVector::add(const SVector &v2, double c2, const FVector &q2) {
	w.detach();
	Rep *r = rep();
	int m = r->size;
	m = max(m, v2.size());
	m = min(m, q2.size());
	if (m > r->size)
		r->resize(m);
	VFloat *d = r->data;
	const VFloat *q = (const VFloat*) q2;
	int npairs = v2.npairs();
	const SVector::Pair *pairs = v2;
	for (int i = 0; i < npairs && pairs->i < m; i++, pairs++)
		d[pairs->i] += pairs->v * q[pairs->i] * c2;
}

void FVector::scale(double c1) {
	w.detach();
	Rep *r = rep();
	VFloat *d = r->data;
	VFloat c = c1;
	int m = r->size;
#if defined(__GNUC__) && defined(__SSE2__) && !defined(__APPLE__)
	typedef float v4sf __attribute__ ((vector_size (16)));
	if (sizeof(VFloat) == 4 && m >= 16)
	{
		VFloat cc[4];
		cc[0] = cc[1] = cc[2] = cc[3] = c;
		v4sf cv = *(v4sf*)cc;
		while (m >= 4)
		{
			v4sf dv = *(v4sf*)d;
			*(v4sf*)d = dv * cv;
			d += 4;
			m -= 4;
		}
	}
#endif
	while (--m >= 0)
		(*d++) *= c;
}

void FVector::combine(double c1d, const FVector &v2, double c2d) {
	w.detach();
	Rep *r = rep();
	int m = max(r->size, v2.size());
	if (m > r->size)
		r->resize(m);
	VFloat *d = r->data;
	const VFloat *s = (const VFloat*) v2;
	VFloat c1 = c1d;
	VFloat c2 = c2d;
#if defined(__GNUC__) && defined(__SSE2__) && !defined(__APPLE__)
	typedef float v4sf __attribute__ ((vector_size (16)));
	if (sizeof(VFloat) == 4 && m >= 16)
	{
		VFloat cc[4];
		cc[0] = cc[1] = cc[2] = cc[3] = c1;
		v4sf c1v = *(v4sf*)cc;
		cc[0] = cc[1] = cc[2] = cc[3] = c2;
		v4sf c2v = *(v4sf*)cc;
		while (m >= 4)
		{
			v4sf sv = *(v4sf*)s;
			v4sf dv = *(v4sf*)d;
			*(v4sf*)d = dv * c1v + sv * c2v;
			d += 4;
			s += 4;
			m -= 4;
		}
	}
#endif
	for (int i = 0; i < m; i++)
		d[i] = d[i] * c1 + s[i] * c2;
}

void FVector::combine(double c1, const SVector &v2, double c2) {
	w.detach();
	Rep *r = rep();
	int m = max(r->size, v2.size());
	if (m > r->size)
		r->resize(m);
	VFloat *d = r->data;
	int npairs = v2.npairs();
	const SVector::Pair *pairs = v2;
	int j = 0;
	for (int i = 0; i < npairs; i++, pairs++) {
		if (c1 != 1)
			for (; j < pairs->i; j++)
				d[j] = d[j] * c1;
		j = pairs->i;
		d[j] = d[j] * c1 + pairs->v * c2;
		j++;
	}
	if (c1 != 1)
		for (; j < pairs->i; j++)
			d[j] = d[j] * c1;
}

std::ostream&
operator<<(std::ostream &f, const FVector &v) {
	std::ios::fmtflags oldflags = f.flags();
	std::streamsize oldprec = f.precision();
	f << std::scientific << std::setprecision(sizeof(VFloat) == 4 ? 7 : 16);
	for (int i = 0; i < v.size(); i++)
		if (v[i] || (i + 1 == v.size())) {
			f << " " << i << ":";
			VFloat x = v[i];
			short ix = (int) x;
			if (x == (VFloat) ix)
				f << ix;
			else
				f << x;
		}
	f << std::endl;
	f.precision(oldprec);
	f.flags(oldflags);
	return f;
}

std::istream&
operator>>(std::istream &f, FVector &v) {
	int sz = 0;
	int msz = 1024;
	v.clear();
	v.resize(msz);
	for (;;) {
		int c = f.get();
		if (!f.good() || (c == '\n' || c == '\r'))
			break;
		if (::isspace(c))
			continue;
		int i;
		f.unget();
		f >> std::skipws >> i >> std::ws;
		if (f.get() != ':') {
			f.unget();
			f.setstate(std::ios::badbit);
			break;
		}
		double x;
		f >> std::skipws >> x;
		if (!f.good())
			break;
		if (i >= sz)
			sz = i + 1;
		if (i >= msz) {
			while (i >= msz)
				msz += msz;
			v.resize(msz);
		}
		v.set(i, x);
	}
	v.resize(sz);
	return f;
}

bool FVector::save(std::ostream &f) const {
	int i = size();
	const VFloat *d = rep()->data;
	f.write((const char*) &i, sizeof(int));
	f.write((const char*) d, sizeof(VFloat) * i);
	return f.good();
}

bool FVector::load(std::istream &f) {
	int i = 0;
	clear();
	f.read((char*) &i, sizeof(int));
	if (i < 0)
		f.setstate(std::ios::badbit);
	if (!f.good())
		return false;
	resize(i);
	VFloat *d = rep()->data;
	f.read((char*) d, sizeof(VFloat) * i);
	return f.good();
}

// ----------------------------------------

void SVector::Rep::resize(int n) {
	if (n != mpairs) {
		Pair *p = new Pair[n + 1];
		int m = min(n, npairs);
		int i = 0;
		for (; i < m; i++)
			p[i] = pairs[i];
		for (; i <= n; i++)
			p[i].i = -1;
		delete[] pairs;
		pairs = p;
		npairs = m;
		mpairs = n;
		size = (m > 0) ? p[m - 1].i + 1 : 0;
	}
}

SVector::Rep *
SVector::Rep::copy() {
	int n = npairs;
	Pair *p = new Pair[n + 1];
	for (int i = 0; i < n; i++)
		p[i] = pairs[i];
	p[n].i = -1;
	Rep *nr = new Rep;
	delete[] nr->pairs;
	nr->pairs = p;
	nr->size = size;
	nr->npairs = nr->mpairs = n;
	return nr;
}

inline double SVector::Rep::qset(int i, double v) {
	assert(i >= size);
	if (npairs >= mpairs) {
		resize(max(16, mpairs + mpairs));
		assert(npairs < mpairs);
	}
	Pair *p = &pairs[npairs++];
	size = i + 1;
	p->i = i;
	p->v = v;
	return v;
}

SVector::SVector() {
	trim();
}

SVector::SVector(const FVector &v) {
	int m = v.size();
	const VFloat *f = v;
	Rep *r = rep();
	r->resize(m);
	for (int i = 0; i < m; i++)
		if (f[i] != 0)
			r->qset(i, f[i]);
	trim();
}


SVector::SVector(const SVector::Rep &r)
{
  Rep *re = rep();
  re->resize(r.npairs);

  for (int i=0; i<r.npairs; i++)
      re->qset(r.pairs[i].i,r.pairs[i].v);
  trim();
}



std::pair<unsigned, double> SVector::extract_component(unsigned aFeatureOrderIndex) {
	const SVector &v = *this;
	const SVector::Rep *r = v.rep();
	const SVector::Pair *pairs = r->pairs;
	int npairs = r->npairs;
	assert(npairs>=(int)aFeatureOrderIndex);
	std::pair<unsigned, double> p = std::make_pair(pairs[aFeatureOrderIndex].i, pairs[aFeatureOrderIndex].v);
	return p;
}

int SVector::sparse_size() {
	const SVector &v = *this;
	const SVector::Rep *r = v.rep();
	int npairs = r->npairs;
	return npairs;
}

namespace {

SVector::Pair *
search(SVector::Pair *pairs, int npairs, int i) {
	int lo = 0;
	int hi = npairs - 1;
	while (lo <= hi) {
		int d = (lo + hi + 1) / 2;
		if (i == pairs[d].i)
			return &pairs[d];
		else if (i < pairs[d].i)
			hi = d - 1;
		else
			lo = d + 1;
	}
	return 0;
}

}

double SVector::get(int i) const {
	const Rep *r = rep();
	if (i >= r->size)
		return 0;
	// binary search
	assert(i >= 0);
	SVector::Pair *pair = search(r->pairs, r->npairs, i);
	if (pair)
		return pair->v;
	return 0;
}

double SVector::set(int i, double v) {
	w.detach();
	Rep *r = rep();
	if (v) {
		if (i >= r->size)
			return r->qset(i, v);
		assert(i >= 0);
		SVector::Pair *p = search(r->pairs, r->npairs, i);
		if (p)
			return p->v = v;
		if (r->npairs >= r->mpairs) {
			r->resize(max(16, r->mpairs + min(r->mpairs, 4096)));
			assert(r->npairs < r->mpairs);
		}
		SVector::Pair *s = r->pairs;
		p = s + r->npairs;
		r->npairs += 1;
		for (; p > s && p[-1].i > i; p--)
			p[0] = p[-1];
		p[0].i = i;
		p[0].v = v;
	} else {
		SVector::Pair *s = r->pairs;
		SVector::Pair *p = search(s, r->npairs, i);
		if (p) {
			r->npairs -= 1;
			for (; p->i >= 0; p++)
				p[0] = p[1];
		}
	}
	return v;
}

void SVector::clear() {
	w.detach();
	rep()->resize(0);
}

void SVector::trim() {
	w.detach();
	Rep *r = rep();
	r->resize(r->npairs);
}

SVector SVector::slice(int fi, int ti) const {
	assert(ti >= 0);
	assert(ti >= fi);
	SVector y;
	for (Pair *p = rep()->pairs; p->i >= 0 && p->i <= ti; p++)
		if (p->i >= fi)
			y.set(p->i, p->v);
	return y;
}

void SVector::normalize() {
	double norm = dot(*this, *this);
	scale(1 / sqrt(norm));
}

void SVector::add(const SVector &v2) {
	operator=(sgd::combine(*this, 1.0, v2, 1.0));
}

void SVector::add(const SVector &v2, double c2) {
	operator=(sgd::combine(*this, 1.0, v2, c2));
}

void SVector::combine(double c1, const SVector &v2, double c2) {
	operator=(sgd::combine(*this, c1, v2, c2));
}

void SVector::scale(double c1) {
	if (c1) {
		w.detach();
		Rep *r = rep();
		Pair *pairs = r->pairs;
		int npairs = r->npairs;
		for (int i = 0; i < npairs; i++)
			pairs[i].v *= c1;
	} else {
		clear();
	}
}

void SVector::binarize() {
	w.detach();
	Rep *r = rep();
	Pair *pairs = r->pairs;
	int npairs = r->npairs;
	for (int i = 0; i < npairs; i++)
		pairs[i].v = 1;
}

std::vector<std::pair<int, double> > SVector::unpack() {
	std::vector<std::pair<int, double> > vec;
	w.detach();
	Rep *r = rep();
	Pair *pairs = r->pairs;
	int npairs = r->npairs;
	for (int i = 0; i < npairs; i++)
		vec.push_back(std::make_pair(pairs[i].i, pairs[i].v));
	return vec;
}

std::ostream&
operator<<(std::ostream &f, const SVector &v) {
	const SVector::Rep *r = v.rep();
	const SVector::Pair *pairs = r->pairs;
	int npairs = r->npairs;
	std::ios::fmtflags oldflags = f.flags();
	std::streamsize oldprec = f.precision();
	f << std::scientific << std::setprecision(sizeof(VFloat) == 4 ? 7 : 16);
	for (int i = 0; i < npairs; i++) {
		f << " " << pairs[i].i << ":";
		VFloat v = pairs[i].v;
		short iv = (int) v;
		if (v == (VFloat) iv)
			f << iv;
		else
			f << v;
	}
	f << std::endl;
	f.precision(oldprec);
	f.flags(oldflags);
	return f;
}

std::istream&
operator>>(std::istream &f, SVector &v) {
	v.clear();
	for (;;) {
		int c = f.get();
		if (!f.good() || (c == '\n' || c == '\r'))
			break;
		if (::isspace(c))
			continue;
		int i;
		f.unget();
		f >> std::skipws >> i >> std::ws;
		if (f.get() != ':') {
			f.unget();
			f.setstate(std::ios::badbit);
			break;
		}
		double x;
		f >> std::skipws >> x;
		if (!f.good())
			break;
		v.set(i, x);
	}
	v.trim();
	return f;
}

bool SVector::save(std::ostream &f) const {
	const Rep *r = rep();
	const Pair *pairs = r->pairs;
	int npairs = r->npairs;
	f.write((const char*) &npairs, sizeof(int));
	f.write((const char*) pairs, sizeof(Pair) * npairs);
	return f.good();
}

bool SVector::load(std::istream &f) {
	clear();
	int npairs = 0;
	f.read((char*) &npairs, sizeof(int));
	if (npairs < 0)
		f.setstate(std::ios::badbit);
	if (!f.good())
		return false;
	rep()->resize(npairs);
	for (int i = 0; i < npairs; i++) {
		Pair pair;
		f.read((char*) &pair, sizeof(Pair));
		if (f.good())
			set(pair.i, pair.v);
	}
	trim();
	return f.good();
}

double dot(const FVector &v1, const FVector &v2) {
	int m = min(v1.size(), v2.size());
	const VFloat *f1 = v1;
	const VFloat *f2 = v2;
	VFloat sum = 0.0;
#if defined(__GNUC__) && defined(__SSE2__) && !defined(__APPLE__)
	typedef float v4sf __attribute__ ((vector_size (16)));
	if (sizeof(VFloat) == 4 && m >= 16)
	{
		VFloat sums[4] = {0,0,0,0};
		v4sf acc = *(v4sf*)sums;
		while (m >= 4)
		{
			v4sf q1 = *(v4sf*)f1;
			v4sf q2 = *(v4sf*)f2;
			acc = acc + q1 * q2;
			m -= 4;
			f1 += 4;
			f2 += 4;
		}
		*(v4sf*)sums = acc;
		sum = sums[0] + sums[1] + sums[2] + sums[3];
	}
#endif
	while (--m >= 0)
		sum += (*f1++) * (*f2++);
	return sum;
}

double dot(const FVector &v1, const SVector &v2) {
	int m = v1.size();
	const VFloat *f = v1;
	const SVector::Pair *p = v2;
	double sum = 0;
	if (p)
		for (; p->i >= 0 && p->i < m; p++)
			sum += (double) p->v * (double) f[p->i];
	return sum;
}

double dot(const SVector &v1, const FVector &v2) {
	int m = v2.size();
	const VFloat *f = v2;
	const SVector::Pair *p = v1;
	double sum = 0;
	if (p)
		for (; p->i >= 0 && p->i < m; p++)
			sum += (double) p->v * (double) f[p->i];
	return sum;
}

double dot(const SVector &v1, const SVector &v2) {
	const SVector::Pair *p1 = v1;
	const SVector::Pair *p2 = v2;
	double sum = 0;
	if (p1 && p2)
		while (p1->i >= 0 && p2->i >= 0) {
			if (p1->i < p2->i)
				p1++;
			else if (p1->i > p2->i)
				p2++;
			else
				sum += (double) (p1++)->v * (double) (p2++)->v;
		}
	return sum;
}

SVector combine(const SVector &v1, double a1, const SVector &v2, double a2) {
	const SVector::Pair *p1 = v1;
	const SVector::Pair *p2 = v2;
	SVector ans;
	SVector::Rep *r = ans.rep();
	r->resize(v1.npairs() + v2.npairs());
	SVector::Pair *p = r->pairs;
	while (p1->i >= 0 && p2->i >= 0) {
		if (p1->i < p2->i) {
			double v = (double) p1->v * a1;
			if (v) {
				p->i = p1->i;
				p->v = v;
				p++;
			}
			p1++;
		} else if (p1->i > p2->i) {
			double v = (double) p2->v * a2;
			if (v) {
				p->i = p2->i;
				p->v = v;
				p++;
			}
			p2++;
		} else {
			double v = (double) p1->v * a1 + (double) p2->v * a2;
			if (v) {
				p->i = p1->i;
				p->v = v;
				p++;
			}
			p1++;
			p2++;
		}
	}
	while (p1->i >= 0) {
		double v = (double) p1->v * a1;
		if (v) {
			p->i = p1->i;
			p->v = v;
			p++;
		}
		p1++;
	}
	while (p2->i >= 0) {
		double v = (double) p2->v * a2;
		if (v) {
			p->i = p2->i;
			p->v = v;
			p++;
		}
		p2++;
	}
	r->npairs = p - r->pairs;
	r->size = (r->npairs > 0) ? p[-1].i + 1 : 0;
	ans.trim();
	return ans;
}

FVector combine(const FVector &v1, double a1, const SVector &v2, double a2) {
	FVector r = v1;
	r.combine(a1, v2, a2);
	return r;
}

FVector combine(const SVector &v1, double a1, const FVector &v2, double a2) {
	FVector r = v2;
	r.combine(a2, v1, a1);
	return r;
}

FVector combine(const FVector &v1, double a1, const FVector &v2, double a2) {
	FVector r = v1;
	r.combine(a1, v2, a2);
	return r;
}

SVector elementwise_product(const SVector &v1, const SVector &v2, const double hardlimit){
	const SVector::Pair *p1 = v1;
	const SVector::Pair *p2 = v2;
	SVector ans;
	SVector::Rep *r = ans.rep();
	r->resize(min(v1.npairs(), v2.npairs()));
	SVector::Pair *p = r->pairs;
	while (p1->i >= 0 && p2->i >= 0) {
		if (p1->i < p2->i) {
			p1++;
		} else if (p1->i > p2->i) {
			p2++;
		} else {
			double v = (double) p1->v * (double) p2->v;
			if (v > hardlimit) v=hardlimit;
			if (v < -hardlimit) v=-hardlimit;
			if (v) {
				p->i = p1->i;
				p->v = v;
				p++;
			}
			p1++;
			p2++;
		}
	}
	r->npairs = p - r->pairs;
	r->size = (r->npairs > 0) ? p[-1].i + 1 : 0;
	ans.trim();
	return ans;
}

} // namespace


/* -------------------------------------------------------------
   Local Variables:
   c++-font-lock-extra-types: ( "\\sw+_t" "[A-Z]\\sw*[a-z]\\sw*" "std::\\sw+")
   End:
   ------------------------------------------------------------- */

