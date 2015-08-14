/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifndef _copyrite_h
#define _copyrite_h

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Has to be a function, so we can get version number from the build system */
const char *GromacsVersion(void);

const char *Program(void);

const char *ShortProgram(void);

/* For both bromacs() and cool_quote() you have to provide a pointer to
 * a string of reasonable length (say 256) and the string length. This
 * is necessary to make the routines threadsafe and avoid allocating
 * a new string each time. The retstring pointer will be the return value.
 */
void
bromacs(char *retstring, int retsize);

/* For cool_quote, the number of the quote used will be returned in cqnum
 * if it is non-NULL.
 */
void
cool_quote(char *retstring, int retsize, int *cqnum);

void
gmx_thanx(FILE *fp);

void
please_cite(FILE *fp, const char *key);
/* Print a message asking to cite something... */

#ifdef __cplusplus
}

namespace gmx
{

class IProgramContext;

/*! \brief
 * Settings for printBinaryInformation().
 *
 * This class is used to specify what printBinaryInformation() prints.
 */
class BinaryInformationSettings
{
    public:
        BinaryInformationSettings();

        //! Print information about build settings.
        BinaryInformationSettings &extendedInfo(bool bEnabled)
        {
            bExtendedInfo_ = bEnabled;
            return *this;
        }
        //! Print copyright and license information.
        BinaryInformationSettings &copyright(bool bEnabled)
        {
            bCopyright_ = bEnabled;
            return *this;
        }
        //! Print a header line with "Generated by" text (for output files).
        BinaryInformationSettings &generatedByHeader(bool bEnabled)
        {
            bGeneratedByHeader_ = bEnabled;
            return *this;
        }
        //! Prefix each line with this string.
        BinaryInformationSettings &linePrefix(const char *prefix)
        {
            prefix_ = prefix;
            return *this;
        }
        //! Suffix each line with this string.
        BinaryInformationSettings &lineSuffix(const char *suffix)
        {
            suffix_ = suffix;
            return *this;
        }

    private:
        bool        bExtendedInfo_;
        bool        bCopyright_;
        bool        bGeneratedByHeader_;
        const char *prefix_;
        const char *suffix_;

        //! Needed to read the members without otherwise unnecessary accessors.
        friend void printBinaryInformation(
            FILE *fp, const IProgramContext &programContext,
            const BinaryInformationSettings &settings);
};

/*! \brief
 * Print basic information about the executable.
 *
 * \param     fp             Where to print the information to.
 * \param[in] programContext Program information object to use.
 */
void printBinaryInformation(FILE                          *fp,
                            const IProgramContext         &programContext);
/*! \brief
 * Print basic information about the executable with custom settings.
 *
 * \param     fp             Where to print the information to.
 * \param[in] programContext Program information object to use.
 * \param[in] settings       Specifies what to print.
 *
 * \see BinaryInformationSettings
 */
void printBinaryInformation(FILE                            *fp,
                            const IProgramContext           &programContext,
                            const BinaryInformationSettings &settings);

} // namespace gmx;

#endif

#endif  /* _copyright_h */