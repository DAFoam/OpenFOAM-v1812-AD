/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "FieldM.H"
#include "FieldFieldReuseFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// CodiPack4OpenFOAM NOTE we need to add Foam:: for Func to prevent
// calling Func from system header
// NOTE: Here Foam::Func(res[i], f[i]) calls function defined in FieldFunctionsM.C
// because res[i] is a Field, e.g., res is a boundaryFields

#define UNARY_FUNCTION(ReturnType, Type, Func)                                 \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    FieldField<Field, ReturnType>& res,                                        \
    const FieldField<Field, Type>& f                                           \
)                                                                              \
{                                                                              \
    forAll(res, i)                                                             \
    {                                                                          \
        Foam::Func(res[i], f[i]);                                              \
    }                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const FieldField<Field, Type>& f                                           \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tRes                                    \
    (                                                                          \
        FieldField<Field, ReturnType>::NewCalculatedType(f)                    \
    );                                                                         \
    Foam::Func(tRes.ref(), f);                                                 \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const tmp<FieldField<Field, Type>>& tf                                     \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tRes(New(tf));                          \
    Foam::Func(tRes.ref(), tf());                                              \
    tf.clear();                                                                \
    return tRes;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_OPERATOR(ReturnType, Type, Op, OpFunc)                           \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    FieldField<Field, ReturnType>& res,                                        \
    const FieldField<Field, Type>& f                                           \
)                                                                              \
{                                                                              \
    forAll(res, i)                                                             \
    {                                                                          \
        Foam::OpFunc(res[i], f[i]);                                            \
    }                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const FieldField<Field, Type>& f                                           \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tRes                                    \
    (                                                                          \
        FieldField<Field, Type>::NewCalculatedType(f)                          \
    );                                                                         \
    Foam::OpFunc(tRes.ref(), f);                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const tmp<FieldField<Field, Type>>& tf                                     \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tRes(New(tf));                          \
    Foam::OpFunc(tRes.ref(), tf());                                            \
    tf.clear();                                                                \
    return tRes;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_FUNCTION(ReturnType, Type1, Type2, Func)                        \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    FieldField<Field, ReturnType>& f,                                          \
    const FieldField<Field, Type1>& f1,                                        \
    const FieldField<Field, Type2>& f2                                         \
)                                                                              \
{                                                                              \
    forAll(f, i)                                                               \
    {                                                                          \
        Foam::Func(f[i], f1[i], f2[i]);                                        \
    }                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const FieldField<Field, Type2>& f2                                         \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tRes                                    \
    (                                                                          \
        FieldField<Field, Type1>::NewCalculatedType(f1)                        \
    );                                                                         \
    Foam::Func(tRes.ref(), f1, f2);                                            \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const tmp<FieldField<Field, Type2>>& tf2                                   \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tRes                                    \
    (                                                                          \
        reuseTmpFieldField<Field, ReturnType, Type2>::New(tf2)                 \
    );                                                                         \
    Foam::Func(tRes.ref(), f1, tf2());                                         \
    tf2.clear();                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const FieldField<Field, Type2>& f2                                         \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tRes                                    \
    (                                                                          \
        reuseTmpFieldField<Field, ReturnType, Type1>::New(tf1)                 \
    );                                                                         \
    Foam::Func(tRes.ref(), tf1(), f2);                                         \
    tf1.clear();                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const tmp<FieldField<Field, Type2>>& tf2                                   \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tRes                                    \
    (                                                                          \
        reuseTmpTmpFieldField<Field, ReturnType, Type1, Type1, Type2>::        \
            New(tf1, tf2)                                                      \
    );                                                                         \
    Foam::Func(tRes.ref(), tf1(), tf2());                                      \
    tf1.clear();                                                               \
    tf2.clear();                                                               \
    return tRes;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    FieldField<Field, ReturnType>& f,                                          \
    const Type1& s,                                                            \
    const FieldField<Field, Type2>& f2                                         \
)                                                                              \
{                                                                              \
    forAll(f, i)                                                               \
    {                                                                          \
        Foam::Func(f[i], s, f2[i]);                                            \
    }                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const Type1& s,                                                            \
    const FieldField<Field, Type2>& f2                                         \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tRes                                    \
    (                                                                          \
        FieldField<Field, Type2>::NewCalculatedType(f2)                        \
    );                                                                         \
    Foam::Func(tRes.ref(), s, f2);                                             \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const Type1& s,                                                            \
    const tmp<FieldField<Field, Type2>>& tf2                                   \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tRes                                    \
    (                                                                          \
        reuseTmpFieldField<Field, ReturnType, Type2>::New(tf2)                 \
    );                                                                         \
    Foam::Func(tRes.ref(), s, tf2());                                          \
    tf2.clear();                                                               \
    return tRes;                                                               \
}


#define BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)                \
                                                                               \
TEMPLATE                                                                       \
void Func                                                                      \
(                                                                              \
    FieldField<Field, ReturnType>& f,                                          \
    const FieldField<Field, Type1>& f1,                                        \
    const Type2& s                                                             \
)                                                                              \
{                                                                              \
    forAll(f, i)                                                               \
    {                                                                          \
        Foam::Func(f[i], f1[i], s);                                            \
    }                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const Type2& s                                                             \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tRes                                    \
    (                                                                          \
        FieldField<Field, Type1>::NewCalculatedType(f1)                        \
    );                                                                         \
    Foam::Func(tRes.ref(), f1, s);                                             \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> Func                                        \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const Type2& s                                                             \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tRes                                    \
    (                                                                          \
        reuseTmpFieldField<Field, ReturnType, Type1>::New(tf1)                 \
    );                                                                         \
    Foam::Func(tRes.ref(), tf1(), s);                                          \
    tf1.clear();                                                               \
    return tRes;                                                               \
}


#define BINARY_TYPE_FUNCTION(ReturnType, Type1, Type2, Func)                   \
    BINARY_TYPE_FUNCTION_SF(ReturnType, Type1, Type2, Func)                    \
    BINARY_TYPE_FUNCTION_FS(ReturnType, Type1, Type2, Func)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)                  \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    FieldField<Field, ReturnType>& f,                                          \
    const FieldField<Field, Type1>& f1,                                        \
    const FieldField<Field, Type2>& f2                                         \
)                                                                              \
{                                                                              \
    forAll(f, i)                                                               \
    {                                                                          \
        Foam::OpFunc(f[i], f1[i], f2[i]);                                      \
    }                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const FieldField<Field, Type2>& f2                                         \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tRes                                    \
    (                                                                          \
        FieldField<Field, ReturnType>::NewCalculatedType(f1)                   \
    );                                                                         \
    Foam::OpFunc(tRes.ref(), f1, f2);                                          \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const tmp<FieldField<Field, Type2>>& tf2                                   \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tRes                                    \
    (                                                                          \
        reuseTmpFieldField<Field, ReturnType, Type2>::New(tf2)                 \
    );                                                                         \
    Foam::OpFunc(tRes.ref(), f1, tf2());                                       \
    tf2.clear();                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const FieldField<Field, Type2>& f2                                         \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tRes                                    \
    (                                                                          \
        reuseTmpFieldField<Field, ReturnType, Type1>::New(tf1)                 \
    );                                                                         \
    Foam::OpFunc(tRes.ref(), tf1(), f2);                                       \
    tf1.clear();                                                               \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const tmp<FieldField<Field, Type2>>& tf2                                   \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tRes                                    \
    (                                                                          \
        reuseTmpTmpFieldField<Field, ReturnType, Type1, Type1, Type2>::        \
            New(tf1, tf2)                                                      \
    );                                                                         \
    Foam::OpFunc(tRes.ref(), tf1(), tf2());                                    \
    tf1.clear();                                                               \
    tf2.clear();                                                               \
    return tRes;                                                               \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)          \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    FieldField<Field, ReturnType>& f,                                          \
    const Type1& s,                                                            \
    const FieldField<Field, Type2>& f2                                         \
)                                                                              \
{                                                                              \
    forAll(f, i)                                                               \
    {                                                                          \
        Foam::OpFunc(f[i], s, f2[i]);                                          \
    }                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const Type1& s,                                                            \
    const FieldField<Field, Type2>& f2                                         \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tRes                                    \
    (                                                                          \
        FieldField<Field, Type2>::NewCalculatedType(f2)                        \
    );                                                                         \
    Foam::OpFunc(tRes.ref(), s, f2);                                           \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const Type1& s,                                                            \
    const tmp<FieldField<Field, Type2>>& tf2                                   \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tRes                                    \
    (                                                                          \
        reuseTmpFieldField<Field, ReturnType, Type2>::New(tf2)                 \
    );                                                                         \
    Foam::OpFunc(tRes.ref(), s, tf2());                                        \
    tf2.clear();                                                               \
    return tRes;                                                               \
}


#define BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)          \
                                                                               \
TEMPLATE                                                                       \
void OpFunc                                                                    \
(                                                                              \
    FieldField<Field, ReturnType>& f,                                          \
    const FieldField<Field, Type1>& f1,                                        \
    const Type2& s                                                             \
)                                                                              \
{                                                                              \
    forAll(f, i)                                                               \
    {                                                                          \
        Foam::OpFunc(f[i], f1[i], s);                                          \
    }                                                                          \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const FieldField<Field, Type1>& f1,                                        \
    const Type2& s                                                             \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tRes                                    \
    (                                                                          \
        FieldField<Field, Type1>::NewCalculatedType(f1)                        \
    );                                                                         \
    Foam::OpFunc(tRes.ref(), f1, s);                                           \
    return tRes;                                                               \
}                                                                              \
                                                                               \
TEMPLATE                                                                       \
tmp<FieldField<Field, ReturnType>> operator Op                                 \
(                                                                              \
    const tmp<FieldField<Field, Type1>>& tf1,                                  \
    const Type2& s                                                             \
)                                                                              \
{                                                                              \
    tmp<FieldField<Field, ReturnType>> tRes                                    \
    (                                                                          \
        reuseTmpFieldField<Field, ReturnType, Type1>::New(tf1)                 \
    );                                                                         \
    Foam::OpFunc(tRes.ref(), tf1(), s);                                        \
    tf1.clear();                                                               \
    return tRes;                                                               \
}


#define BINARY_TYPE_OPERATOR(ReturnType, Type1, Type2, Op, OpFunc)             \
    BINARY_TYPE_OPERATOR_SF(ReturnType, Type1, Type2, Op, OpFunc)              \
    BINARY_TYPE_OPERATOR_FS(ReturnType, Type1, Type2, Op, OpFunc)

// ************************************************************************* //
