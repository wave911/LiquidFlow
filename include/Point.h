/*
 * Point.h
 *
 *  Created on: Apr 27, 2013
 *      Author: ubuntu
 */

#ifndef POINT_H_
#define POINT_H_

#include "common_headers/real_type.h"

class CPoint3D
{
	public:
		real_t m_x,
			   m_y,
			   m_z;

	public:
		CPoint3D() {
			m_x = 0;
			m_y = 0;
			m_z = 0;
		};

		CPoint3D(real_t x, real_t y, real_t z)
		{
			m_x = x;
			m_y = y;
			m_z = z;
		};

		CPoint3D(const CPoint3D &pt)
		{
			*this = pt;
		};

		virtual ~CPoint3D() {};

		bool operator== (const CPoint3D &pt)
		{
			if ((this->m_x == pt.m_x) && (this->m_y == pt.m_y) && (this->m_z == pt.m_z))
				return true;
			else
				return false;
		};

		CPoint3D &operator= (const CPoint3D &pt)
		{
			this->m_x = pt.m_x;
			this->m_y = pt.m_y;
			this->m_z = pt.m_z;

			return *this;
		}
};

class CPointProperties
{
	public:
		int m_num;
		bool m_isborder;
		CPoint3D m_norm;
	public:
		CPointProperties()
		{
			m_num = 0;
			m_isborder = false;
			m_norm.m_x = 0; m_norm.m_y = 0; m_norm.m_z =0;
		};

		CPointProperties(int num, bool isborder)
		{
			m_num = num;
			m_isborder = isborder;
		};

		CPointProperties(int num, bool isborder, CPoint3D norm)
		{
			m_num = num;
			m_isborder = isborder;
			m_norm = norm;
		};

		~CPointProperties() {};
};


#endif /* POINT_H_ */
